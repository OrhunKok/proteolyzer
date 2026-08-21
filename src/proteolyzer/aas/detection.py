"""Detection of candidate amino acid substitutions and PTMs.

:class:`Detection` is a thin dispatcher: it reads the ``Utils.Workflow``
parameter and delegates to the matching subclass (currently only
:class:`MaxQuant`).
"""

import pickle
import re
import shutil

import ahocorasick
import numpy as np
import pandas as pd

from proteolyzer import reference
from proteolyzer.core.io import read_frame, write_frame

from . import utils
from .base import Stage
from .config import Config

CONFIG = Config()


class Detection(Stage):
    """Dispatches to the search-engine specific detection workflow."""

    def __init__(self, params, queue=None):
        super().__init__(params, queue)

        self.subs_ref = utils.aa_subs_ref()
        self.mods_dict = utils.gen_mod_dict()
        self.protease = self.params["Detection"]["Protease"]
        self.prot_fasta = self.params["Detection"]["Protein FASTA"]
        self.translated_frames = self.params["Translation"]["Translated Frames Folder"]

        if type(self) is Detection:
            self._initialize_workflow()

    def _initialize_workflow(self):
        """Instantiate the subclass named by the configured workflow."""
        workflows = {cls.__name__: cls for cls in type(self).__subclasses__()}
        workflow = workflows.get(self.workflow)

        if not workflow:
            raise NotImplementedError(
                f"Detection '{self.workflow}' class is not found. "
                f"Available workflows: {sorted(workflows)}"
            )

        self.detection_workflow = workflow(self._raw_params, self.queue)

    def run(self):
        """Run the detection workflow."""
        self.detection_workflow.run()

    def process_sample(self, sample):
        """Delegate sample processing to the detection workflow."""
        self.detection_workflow.process_sample(sample)


class MaxQuant(Detection):
    def __init__(self, params, queue=None):
        super().__init__(params, queue)

        self.mq_params = self.params["Detection"]
        # "Detection PEP" is applied by Preprocessor.MaxQuant, which is where
        # dependent peptides are filtered.
        self.aa_sub_ppm = float(self.mq_params["AA Substitution ppm"])
        self.ptm_ppm = float(self.mq_params["PTM ppm"])
        self.pos_prob = float(self.mq_params["Positional Probability Threshold"])
        self.cn_term_prob = float(self.mq_params["C/n-term Modification Threshold"])

    def run(self):
        """Run the MaxQuant-specific detection workflow over every sample."""
        self.record_run()
        for sample in self.samples:
            self.process_sample(sample)

    def process_sample(self, sample):
        """Process each sample with MaxQuant-specific logic."""
        sample_dir = self._locate_sample_dir(sample)
        if not sample_dir:
            self.queue.put(("stderr", f"Sample {sample} not found in data directory."))
            return

        self.queue.put(("stdout", f"Processing sample: {sample}"))

        all_peps = pd.read_parquet(
            sample_dir / "allPeptides.parquet", engine="fastparquet"
        )
        peptides = pd.read_parquet(
            sample_dir / "peptides.parquet", engine="fastparquet"
        )

        n_term_peps = peptides.loc[peptides["Terminus"] == "N", "Sequence"].to_numpy()
        c_term_peps = peptides.loc[peptides["Terminus"] == "C", "Sequence"].to_numpy()

        sample_df = self.find_potential_aas(all_peps, n_term_peps, c_term_peps)
        utils.ptm_mtp_output(sample_df, sample, self.output_dir)
        self._validate_mtp(sample, sample_dir / "evidence.parquet")

    def _validate_mtp(self, sample, evidence_path):
        """Validate MTP and output filtered MTP."""
        mtp = read_frame(self.output_dir / "MTP" / f"{sample}_MTP")

        evidence = pd.read_parquet(evidence_path, engine="fastparquet")
        mtp = self.validate_mtp(evidence, mtp)

        write_frame(mtp, self.output_dir / "MTP" / f"{sample}_MTP_Filtered_Stage_1")

        self.write_fasta(mtp, sample)

    def find_potential_aas(self, dp_df, n_term_peps, c_term_peps):
        """Identify potential AA substitutions and other modifications."""
        dp_df = self._prepare_peptide_data(dp_df, n_term_peps, c_term_peps)
        dp_df = self._apply_modifications(dp_df)
        dp_df = self._find_homologous_peptides(dp_df)
        return dp_df

    def _prepare_peptide_data(self, dp_df, n_term_peps, c_term_peps):
        """Prepare the peptide data for further processing."""
        dp_df["Leading.Razor.DP.Protein"] = dp_df["DP Proteins"].str.split(";").str[0]
        dp_df["Con.DP.Protein"] = dp_df["DP Proteins"].str.contains("CON__", na=False)

        dp_df["DP Positions"] = dp_df.apply(
            lambda x: [i - 1 for i, a in enumerate(x["DP Probabilities"]) if a == "("],
            axis=1,
        )
        (
            dp_df["count candidate residues per peptide"],
            dp_df["DP candidate residues"],
            dp_df["DP positional probabilities"],
        ) = zip(
            *dp_df.apply(
                lambda x: self.refine_localization_probabilities(x["DP Probabilities"]),
                axis=1,
            ),
            strict=True,
        )
        dp_df["Protein N-term"] = dp_df["DP Base Sequence"].isin(n_term_peps)
        dp_df["Protein C-term"] = dp_df["DP Base Sequence"].isin(c_term_peps)
        dp_df["Peptide N-term"] = dp_df["DP Probabilities"].apply(self.pep_nterm)
        dp_df["Peptide C-term"] = dp_df["DP Probabilities"].apply(self.pep_cterm)
        dp_df = dp_df.explode(
            ["DP candidate residues", "DP positional probabilities", "DP Positions"]
        )
        self.queue.put(("progress", (1, 4)))
        return dp_df

    def _apply_modifications(self, dp_df):
        """Apply AA substitutions and PTM modifications to peptides."""
        dp_df[
            [
                "aa subs",
                "aa subs positional probability",
                "aa subs mass error (ppm)",
                "destination aa",
            ]
        ] = None
        dp_df = dp_df.apply(lambda x: self.get_aa_subs(x), axis=1)

        dp_df[["mistranslated sequence", "mistranslated aas positions"]] = None
        dp_df = dp_df.apply(
            lambda x: self.get_mistranslated_seq(x) if x["aa subs"] else x, axis=1
        )

        dp_df[
            [
                "PTM",
                "PTM site",
                "PTM positional probability",
                "PTM mass error [observed-expected] (ppm)",
            ]
        ] = None
        dp_df = dp_df.apply(lambda x: self.find_PTMs(x), axis=1)
        self.queue.put(("progress", (2, 4)))
        return dp_df

    def refine_localization_probabilities(self, modified_seq):
        modified_sites = [
            modified_seq[m.start() - 1] for m in re.finditer(r"\(", modified_seq)
        ]
        weights = [float(i) for i in re.findall(r"\(([^)]+)\)", modified_seq)]
        return len(modified_sites), modified_sites, weights

    def pep_cterm(self, modified_sequence):
        if modified_sequence[-1] == ")":
            prob = float(modified_sequence[:-1].split("(")[-1])
            return prob >= self.cn_term_prob
        return False

    def pep_nterm(self, modified_sequence):
        if modified_sequence[1] == "(":
            prob = float(modified_sequence[2:].split(")")[0])
            return prob >= self.cn_term_prob
        return False

    def _find_homologous_peptides(self, dp_df):
        """Find homologous peptides based on translated frames."""
        dp_df = self.find_homologous_peptide(dp_df)
        self.queue.put(("progress", (3, 4)))
        return dp_df

    def get_aa_subs(self, row):
        """Get the amino acid substitutions for a given row."""
        candidate = row["DP candidate residues"]
        pos_probs = row["DP positional probabilities"]
        DP_deltam = row["DP Mass Difference"]
        mtol = row["m/z"] * (self.aa_sub_ppm / 1e6)

        res_dict = self.subs_ref.get(candidate, {})
        for s in res_dict:
            delta_m = res_dict[s]
            if (
                ((DP_deltam > 0 and delta_m > 0) or (DP_deltam < 0 and delta_m < 0))
                and (DP_deltam > delta_m - mtol)
                and (DP_deltam < delta_m + mtol)
                and (self.pos_prob < pos_probs)
            ):
                row["aa subs"] = s
                row["aa subs positional probability"] = pos_probs
                row["aa subs mass error (ppm)"] = abs(DP_deltam - delta_m)
                row["destination aa"] = s[-1]
        return row

    def get_mistranslated_seq(self, row):
        """Generate the mistranslated sequence based on substitutions."""
        bp = row["DP Base Sequence"]
        seq = row["DP Probabilities"]
        sub = row["destination aa"]
        pos = row["DP Positions"]
        parsed_seq = np.array(list(seq))
        mtp = parsed_seq.copy()
        mtp[pos] = sub
        mtp = "".join(mtp)
        mtp = re.sub(r"\([^)]*\)", "", mtp)

        position = int([i for i, x in enumerate(bp) if mtp[i] != x][0])
        row["mistranslated sequence"] = mtp
        row["mistranslated aas positions"] = position
        return row

    def find_PTMs(self, row):
        """Find potential PTMs for a given peptide."""
        res = row["DP candidate residues"]
        prob = row["DP positional probabilities"]
        DP_deltam = row["DP Mass Difference"]
        mtol = row["m/z"] * (self.ptm_ppm / 1e6)

        if row["Peptide N-term"]:
            res = "N-term"
        elif row["Peptide C-term"]:
            res = "C-term"

        filtered_mods = self.mods_dict.get(res, [])
        for mod in filtered_mods:
            modification, pos, delta_m = mod
            if (DP_deltam > delta_m - mtol) and (DP_deltam < delta_m + mtol):
                term_filter = (
                    pos == "Anywhere"
                    or (pos == "Protein N-term" and row["Protein N-term"])
                    or (pos == "Any N-term" and row["Peptide N-term"])
                    or (pos == "Protein C-term" and row["Protein C-term"])
                    or (pos == "Any C-term" and row["Peptide C-term"])
                )
                if term_filter:
                    row["PTM"] = modification
                    row["PTM site"] = res
                    row["PTM positional probability"] = prob
                    row["PTM mass error [observed-expected] (ppm)"] = (
                        DP_deltam - delta_m
                    )
        return row

    def aho_corasick_search(self, translated_genome, patterns):
        """Search for patterns in the translated genome using Aho-Corasick."""
        automaton = ahocorasick.Automaton()
        results = {pattern: [] for pattern in patterns}
        for pattern in patterns:
            automaton.add_word(pattern, pattern)
        automaton.make_automaton()

        for end_index, pattern in automaton.iter(translated_genome):
            results[pattern].append(end_index - len(pattern) + 1)

        return {p: (bool(v), v) for p, v in results.items()}

    def aho_corasick_output_organize(self, result):
        """Organize the output from the Aho-Corasick search."""
        return [k[1:] for k, v in result.items() if v[0]]

    def find_homologous_peptide(self, sample_df):
        """Find homologous peptides using the Aho-Corasick algorithm."""
        # A new list, not `+=`: CLEAVAGE_SITES is shared class-level state.
        # The stop codon is a valid preceding residue here too.
        cleavage_sites = [*reference.protease(self.protease).cleavage_sites, "*"]
        all_mtps = sample_df["mistranslated sequence"].dropna().unique()
        all_mtps = [prefix + seq for seq in all_mtps for prefix in cleavage_sites]

        for frame in range(1, 7):
            with open(
                self.translated_frames / f"frame_{frame}_il_ambigous.p", "rb"
            ) as wf:
                w_aa = pickle.load(wf)
            with open(self.translated_frames / f"frame_{frame}.p", "rb") as sf:
                s_aa = pickle.load(sf)
            w_aa_out = self.aho_corasick_search(w_aa, all_mtps)
            s_aa_out = self.aho_corasick_search(s_aa, all_mtps)
            matched = self.aho_corasick_output_organize(
                w_aa_out
            ) + self.aho_corasick_output_organize(s_aa_out)
            col = f"{frame}-frame genome substring"
            sample_df[col] = False
            sample_df.loc[sample_df["mistranslated sequence"].isin(matched), col] = True
        return sample_df

    def evidence_ppm(self, evidence):
        """Calculate the average and standard deviation of mass error (ppm)."""
        me = [x for x in evidence["Mass error [ppm]"] if x > -1000]
        return np.mean(me), np.std(me)

    def posterior_aasub_prob(self, row, mean, std):
        """Calculate posterior probability for amino acid substitution."""
        pp = row["aa subs positional probability"]
        merr = row["aa subs mass error (ppm)"]
        coeff = 1 / (std * np.sqrt(2 * np.pi))
        exponent = np.exp(-0.5 * ((merr - mean) / std) ** 2)
        pdf = coeff * exponent

        return pdf * pp

    def q_val_calc(self, mtp):
        """Calculate q-values for peptides based on posterior probabilities."""
        post_prob = mtp["Posterior subs probability"].to_numpy()
        pval = 1 - post_prob
        ranked_pval = np.sort(pval)
        cumsum = np.cumsum(ranked_pval)
        ranked_qval = np.array([x / (i + 1) for i, x in enumerate(cumsum)])
        qval = ranked_qval[np.argsort(np.argsort(pval))]
        return qval

    def gen_metrics(self, mtp):
        """Generate evaluation metrics."""
        qval = mtp["q-value"].to_numpy()
        thresh = np.max(qval)
        n_thresh = len([x for x in qval if x <= thresh])
        TP_thresh = np.floor((1 - thresh) * n_thresh)
        ref_model = np.logspace(-20, -1, num=100, base=10)
        metric_arr = np.zeros((len(ref_model), 8))
        for i, qt in enumerate(ref_model):
            TP = len([x for x in qval if x <= qt]) * (1 - qt)
            FP = len([x for x in qval if x <= qt]) * qt
            FN = TP_thresh - TP
            TN = len(qval) - TP - FN - FP
            P = TP / (TP + FP) if (TP + FP) > 0 else 0
            R = TP / (TP + FN) if (TP + FN) > 0 else 0
            F_score = (2 * P * R) / (P + R) if (P + R) > 0 else 0
            metric_arr[i, :] = qt, TP, FP, FN, TN, P, R, F_score
        return pd.DataFrame(
            metric_arr,
            columns=[
                "q_threshold",
                "TP",
                "FP",
                "FN",
                "TN",
                "Precision",
                "Recall",
                "F_score",
            ],
        )

    def mtp_filter(self, mtp, metric_df):
        """Filter MTP based on F-score."""
        max_F_idx = metric_df["F_score"].idxmax()
        q_thresh = metric_df.loc[max_F_idx, "q_threshold"]
        filtered_mtp = mtp[mtp["q-value"] <= q_thresh].reset_index(drop=True)
        return filtered_mtp

    def validate_mtp(self, evidence, mtp):
        """Validate MTP using posterior probability and q-value filtering."""
        group_cols = ["Raw file", "DP Base Sequence", "Charge", "DP PEP"]

        # A sample can legitimately end up with no candidates (every mass shift
        # explained by a PTM or by the genome). The q-value model has nothing to
        # rank in that case, so return early rather than reduce empty arrays.
        if mtp.empty:
            self.queue.put(("stdout", "No MTP candidates to validate."))
            return mtp

        mean, std = self.evidence_ppm(evidence)
        mtp["Posterior subs probability"] = mtp.apply(
            lambda x: self.posterior_aasub_prob(x, mean, std), axis=1
        )
        pp_sum = (
            mtp.groupby(group_cols)["Posterior subs probability"].sum().reset_index()
        )
        pp_sum = pp_sum.rename(
            {"Posterior subs probability": "Posterior subs probability sum"}, axis=1
        )
        mtp = mtp.merge(pp_sum, on=group_cols, how="left")
        mtp["Posterior subs probability"] = (
            mtp["Posterior subs probability"]
            / mtp["Posterior subs probability sum"]
            * (1 - mtp["DP PEP"])
        )
        mtp["q-value"] = self.q_val_calc(mtp)
        metrics = self.gen_metrics(mtp)
        return self.mtp_filter(mtp, metrics)

    def write_fasta(self, filtered_mtp, sample_name):
        """Write filtered MTP results to a FASTA file."""

        output_fasta_path = self.output_dir / f"{sample_name}_validation.fasta"
        shutil.copy(self.prot_fasta, output_fasta_path)

        fasta_df = filtered_mtp.loc[
            :,
            [
                "DP Base Sequence",
                "mistranslated sequence",
                "destination aa",
                "mistranslated aas positions",
                "aa subs",
                "Leading.Razor.DP.Protein",
            ],
        ].copy()
        fasta_df["aa subs"] = fasta_df["aa subs"].str.replace(" to ", ":")
        fasta_df["mistranslated aas positions"] = fasta_df[
            "mistranslated aas positions"
        ].astype(int)

        # J (Xle) is ambiguous: emit one entry for isoleucine and one for leucine.
        rows_j = fasta_df.loc[fasta_df["destination aa"] == "J"]
        resolved_j = []
        for aa in ("I", "L"):
            df = rows_j.assign(**{"destination aa": aa})
            df["aa subs"] = df["aa subs"].str[:-1] + aa
            # Built by comprehension rather than DataFrame.apply(axis=1): apply
            # on an empty frame returns a frame, not a series, which cannot be
            # assigned back to a single column.
            df["mistranslated sequence"] = [
                sequence[:position] + aa + sequence[position + 1 :]
                for sequence, position in zip(
                    rows_j["mistranslated sequence"],
                    rows_j["mistranslated aas positions"],
                    strict=True,
                )
            ]
            resolved_j.append(df)

        fasta_df = pd.concat(
            [fasta_df.loc[fasta_df["destination aa"] != "J"], *resolved_j],
            ignore_index=True,
        )

        fasta_df["Header"] = (
            ">MTP|("
            + fasta_df["DP Base Sequence"]
            + ")("
            + fasta_df["mistranslated aas positions"].astype(str)
            + ")("
            + fasta_df["aa subs"]
            + ")("
            + fasta_df["Leading.Razor.DP.Protein"]
            + ")"
        )

        fasta_df = fasta_df.drop_duplicates(subset=["Header"], ignore_index=True)
        write_frame(fasta_df, self.output_dir / "MTP" / f"{sample_name}_FASTA")

        with open(output_fasta_path, "a") as f:
            for header, seq in zip(
                fasta_df["Header"], fasta_df["mistranslated sequence"], strict=True
            ):
                f.write(f"{header}\n{seq}\n")

        self.queue.put(("stdout", "Validation Fasta Written."))
