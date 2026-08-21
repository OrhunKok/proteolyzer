"""Detection of candidate substitutions (SAAP) and their alternatives (ALT).

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
        self.alt_ppm = float(self.mq_params["ALT ppm"])
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
        utils.saap_alt_output(sample_df, sample, self.output_dir)
        self._validate_saap(sample, sample_dir / "evidence.parquet")

    def _validate_saap(self, sample, evidence_path):
        """Validate SAAPs and write the filtered set."""
        saap = read_frame(self.output_dir / "SAAP" / f"{sample}_SAAP")

        evidence = pd.read_parquet(evidence_path, engine="fastparquet")
        saap = self.validate_saap(evidence, saap)

        write_frame(saap, self.output_dir / "SAAP" / f"{sample}_SAAP_Filtered_Stage_1")

        self.write_fasta(saap, sample)

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
        """Assign substitutions, then alternative explanations, to peptides."""
        dp_df[
            [
                "aa subs",
                "aa subs positional probability",
                "aa subs mass error (ppm)",
                "destination aa",
            ]
        ] = None
        dp_df = dp_df.apply(lambda x: self.get_aa_subs(x), axis=1)

        dp_df[["SAAP sequence", "SAAP position"]] = None
        dp_df = dp_df.apply(
            lambda x: self.get_saap_sequence(x) if x["aa subs"] else x, axis=1
        )

        dp_df[
            [
                "ALT",
                "ALT site",
                "ALT positional probability",
                "ALT mass error [observed-expected] (ppm)",
            ]
        ] = None
        dp_df = dp_df.apply(lambda x: self.find_alt(x), axis=1)
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

    def get_saap_sequence(self, row):
        """Generate the SAAP sequence implied by the substitution."""
        base = row["DP Base Sequence"]
        seq = row["DP Probabilities"]
        sub = row["destination aa"]
        pos = row["DP Positions"]
        parsed_seq = np.array(list(seq))
        saap = parsed_seq.copy()
        saap[pos] = sub
        saap = "".join(saap)
        saap = re.sub(r"\([^)]*\)", "", saap)

        position = int([i for i, x in enumerate(base) if saap[i] != x][0])
        row["SAAP sequence"] = saap
        row["SAAP position"] = position
        return row

    def find_alt(self, row):
        """Find an alternative explanation (a known modification) for a peptide."""
        res = row["DP candidate residues"]
        prob = row["DP positional probabilities"]
        DP_deltam = row["DP Mass Difference"]
        mtol = row["m/z"] * (self.alt_ppm / 1e6)

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
                    row["ALT"] = modification
                    row["ALT site"] = res
                    row["ALT positional probability"] = prob
                    row["ALT mass error [observed-expected] (ppm)"] = (
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
        all_saaps = sample_df["SAAP sequence"].dropna().unique()
        all_saaps = [prefix + seq for seq in all_saaps for prefix in cleavage_sites]

        for frame in range(1, 7):
            with open(
                self.translated_frames / f"frame_{frame}_il_ambigous.p", "rb"
            ) as wf:
                w_aa = pickle.load(wf)
            with open(self.translated_frames / f"frame_{frame}.p", "rb") as sf:
                s_aa = pickle.load(sf)
            w_aa_out = self.aho_corasick_search(w_aa, all_saaps)
            s_aa_out = self.aho_corasick_search(s_aa, all_saaps)
            matched = self.aho_corasick_output_organize(
                w_aa_out
            ) + self.aho_corasick_output_organize(s_aa_out)
            col = f"{frame}-frame genome substring"
            sample_df[col] = False
            sample_df.loc[sample_df["SAAP sequence"].isin(matched), col] = True
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

    def q_val_calc(self, saap):
        """Calculate q-values for peptides based on posterior probabilities."""
        post_prob = saap["Posterior subs probability"].to_numpy()
        pval = 1 - post_prob
        ranked_pval = np.sort(pval)
        cumsum = np.cumsum(ranked_pval)
        ranked_qval = np.array([x / (i + 1) for i, x in enumerate(cumsum)])
        qval = ranked_qval[np.argsort(np.argsort(pval))]
        return qval

    def gen_metrics(self, saap):
        """Generate evaluation metrics."""
        qval = saap["q-value"].to_numpy()
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

    def saap_filter(self, saap, metric_df):
        """Filter SAAPs based on F-score."""
        max_F_idx = metric_df["F_score"].idxmax()
        q_thresh = metric_df.loc[max_F_idx, "q_threshold"]
        filtered_saap = saap[saap["q-value"] <= q_thresh].reset_index(drop=True)
        return filtered_saap

    def validate_saap(self, evidence, saap):
        """Validate SAAPs using posterior probability and q-value filtering."""
        group_cols = ["Raw file", "DP Base Sequence", "Charge", "DP PEP"]

        # A sample can legitimately end up with no candidates (every mass shift
        # explained by a modification or by the genome). The q-value model has
        # nothing to rank, so return early rather than reduce empty arrays.
        if saap.empty:
            self.queue.put(("stdout", "No SAAP candidates to validate."))
            return saap

        mean, std = self.evidence_ppm(evidence)
        saap["Posterior subs probability"] = saap.apply(
            lambda x: self.posterior_aasub_prob(x, mean, std), axis=1
        )
        pp_sum = (
            saap.groupby(group_cols)["Posterior subs probability"].sum().reset_index()
        )
        pp_sum = pp_sum.rename(
            {"Posterior subs probability": "Posterior subs probability sum"}, axis=1
        )
        saap = saap.merge(pp_sum, on=group_cols, how="left")
        saap["Posterior subs probability"] = (
            saap["Posterior subs probability"]
            / saap["Posterior subs probability sum"]
            * (1 - saap["DP PEP"])
        )
        saap["q-value"] = self.q_val_calc(saap)
        metrics = self.gen_metrics(saap)
        return self.saap_filter(saap, metrics)

    def write_fasta(self, filtered_saap, sample_name):
        """Write the filtered SAAPs to a FASTA file."""

        output_fasta_path = self.output_dir / f"{sample_name}_validation.fasta"
        shutil.copy(self.prot_fasta, output_fasta_path)

        fasta_df = filtered_saap.loc[
            :,
            [
                "DP Base Sequence",
                "SAAP sequence",
                "destination aa",
                "SAAP position",
                "aa subs",
                "Leading.Razor.DP.Protein",
            ],
        ].copy()
        fasta_df["aa subs"] = fasta_df["aa subs"].str.replace(" to ", ":")
        fasta_df["SAAP position"] = fasta_df["SAAP position"].astype(int)

        # J (Xle) is ambiguous: emit one entry for isoleucine and one for leucine.
        rows_j = fasta_df.loc[fasta_df["destination aa"] == "J"]
        resolved_j = []
        for aa in ("I", "L"):
            df = rows_j.assign(**{"destination aa": aa})
            df["aa subs"] = df["aa subs"].str[:-1] + aa
            # Built by comprehension rather than DataFrame.apply(axis=1): apply
            # on an empty frame returns a frame, not a series, which cannot be
            # assigned back to a single column.
            df["SAAP sequence"] = [
                sequence[:position] + aa + sequence[position + 1 :]
                for sequence, position in zip(
                    rows_j["SAAP sequence"],
                    rows_j["SAAP position"],
                    strict=True,
                )
            ]
            resolved_j.append(df)

        fasta_df = pd.concat(
            [fasta_df.loc[fasta_df["destination aa"] != "J"], *resolved_j],
            ignore_index=True,
        )

        fasta_df["Header"] = (
            ">SAAP|("
            + fasta_df["DP Base Sequence"]
            + ")("
            + fasta_df["SAAP position"].astype(str)
            + ")("
            + fasta_df["aa subs"]
            + ")("
            + fasta_df["Leading.Razor.DP.Protein"]
            + ")"
        )

        fasta_df = fasta_df.drop_duplicates(subset=["Header"], ignore_index=True)
        write_frame(fasta_df, self.output_dir / "SAAP" / f"{sample_name}_FASTA")

        with open(output_fasta_path, "a") as f:
            for header, seq in zip(
                fasta_df["Header"], fasta_df["SAAP sequence"], strict=True
            ):
                f.write(f"{header}\n{seq}\n")

        self.queue.put(("stdout", "Validation Fasta Written."))
