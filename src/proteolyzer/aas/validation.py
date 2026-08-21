"""Fragment-level validation of the substitutions proposed by detection."""

import numpy as np
import pandas as pd

from proteolyzer.core.io import frame_exists, read_frame, write_frame

from .base import Stage


class Validation(Stage):
    """Keeps only candidates supported by fragment ions spanning the substitution."""

    def __init__(self, params, queue=None):
        super().__init__(params, queue)

        # Validation thresholds
        self.validation_pep = self.params["Validation"]["Validation PEP"]
        self.pif = self.params["Validation"]["PIF"]
        self.frag_evidence = self.params["Validation"]["Fragment Evidence"]

    def process_sample(self, sample):
        sample_dir = self._locate_sample_dir(sample, suffix="_val")
        if not sample_dir:
            self.queue.put(("stderr", f"{sample} not found in data directory."))
            return

        self.queue.put(("stdout", f"Processing sample: {sample}"))

        evidence_path = sample_dir / "evidence.parquet"
        msms_path = sample_dir / "msms.parquet"
        saap_path = self.output_dir / "SAAP" / f"{sample}_SAAP_Filtered_Stage_1"

        if not (
            evidence_path.exists() and msms_path.exists() and frame_exists(saap_path)
        ):
            self.queue.put(("stderr", f"Missing data files for sample {sample}"))
            return

        evidence = pd.read_parquet(evidence_path, engine="fastparquet")
        msms = pd.read_parquet(msms_path, engine="fastparquet")
        saap = read_frame(saap_path)

        validated_saap, val_evidence = self.saap_validate(evidence, msms, saap)

        write_frame(
            validated_saap, self.output_dir / "SAAP" / f"{sample}_SAAP_Filtered_Stage_2"
        )
        write_frame(
            val_evidence,
            self.output_dir / "SAAP" / f"{sample}_Val_Evidence_Filtered_Stage_2",
        )

        self.queue.put(("stdout", f"Validation complete for sample: {sample}"))

    def frags_containing_aas(self, row: pd.Series):
        """Fragment indices whose ions span the substituted residue."""
        seq = row["SAAP sequence"]
        pos = row["SAAP position"]

        b_ions_pos = pos + 1
        y_ions_pos = len(seq) - pos

        b_ions_aas = np.arange(b_ions_pos, len(seq), dtype=int)
        y_ions_aas = np.arange(y_ions_pos, len(seq), dtype=int)

        return b_ions_aas, y_ions_aas

    def frag_count(self, row: pd.Series, frag_ev_merge: pd.DataFrame):
        """Count matched b/y fragments supporting the substitution in `row`.

        Kept for use on a single candidate. :meth:`fragment_evidence` computes
        the same numbers for a whole frame without rescanning the fragments
        once per candidate.
        """
        b_ions = row["b_ions_aas"]
        y_ions = row["y_ions_aas"]
        seq = row["SAAP sequence"]
        scan = row["MS/MS scan number"]
        charge = row["Charge"]

        matches = frag_ev_merge[
            (frag_ev_merge["Raw file"] == row["Raw file"])
            & (frag_ev_merge["Sequence"] == seq)
            & (frag_ev_merge["MS/MS scan number"].isin(scan))
            & (frag_ev_merge["Charge"] == charge)
        ]

        b_matches = matches[
            (matches["Frag.Type"] == "b") & (matches["Frag.Number"].isin(b_ions))
        ].shape[0]
        y_matches = matches[
            (matches["Frag.Type"] == "y") & (matches["Frag.Number"].isin(y_ions))
        ].shape[0]

        return b_matches + y_matches

    #: Columns identifying one observed fragment.
    _FRAGMENT_KEYS = (
        "Raw file",
        "Sequence",
        "MS/MS scan number",
        "Charge",
        "Frag.Type",
        "Frag.Number",
    )

    def fragment_evidence(
        self, saap: pd.DataFrame, frag_ev_merge: pd.DataFrame
    ) -> pd.Series:
        """Matched b/y fragment count per candidate, for the whole frame.

        Equivalent to :meth:`frag_count` per row, but by joining the fragments
        a candidate *would* produce against the ones observed. Filtering per
        candidate instead costs one pass over the fragments each, which is
        quadratic in the size of the run.
        """
        wanted = self._expected_fragments(saap)
        if wanted.empty:
            return pd.Series(0, index=saap.index, dtype="int64")

        observed = frag_ev_merge.loc[:, list(self._FRAGMENT_KEYS)].copy()
        for frame in (wanted, observed):
            for key in ("MS/MS scan number", "Charge", "Frag.Number"):
                frame[key] = pd.to_numeric(frame[key], errors="coerce").astype("int64")

        counts = (
            wanted.merge(observed, on=list(self._FRAGMENT_KEYS), how="inner")
            .groupby("_candidate")
            .size()
        )
        return counts.reindex(saap.index, fill_value=0).astype("int64")

    def _expected_fragments(self, saap: pd.DataFrame) -> pd.DataFrame:
        """One row per (candidate, scan, ion) the candidate should produce."""
        keys = pd.DataFrame(
            {
                "_candidate": saap.index,
                "Raw file": saap["Raw file"].to_numpy(),
                "Sequence": saap["SAAP sequence"].to_numpy(),
                "Charge": saap["Charge"].to_numpy(),
                "MS/MS scan number": saap["MS/MS scan number"].to_numpy(),
                "b_ions_aas": saap["b_ions_aas"].to_numpy(),
                "y_ions_aas": saap["y_ions_aas"].to_numpy(),
            }
        ).explode("MS/MS scan number")

        per_type = []
        for frag_type, ion_col in (("b", "b_ions_aas"), ("y", "y_ions_aas")):
            expanded = keys.explode(ion_col)
            expanded = expanded[expanded[ion_col].notna()]
            per_type.append(
                expanded.assign(**{"Frag.Type": frag_type}).rename(
                    columns={ion_col: "Frag.Number"}
                )[["_candidate", *self._FRAGMENT_KEYS]]
            )

        return pd.concat(per_type, ignore_index=True).dropna(
            subset=["MS/MS scan number", "Frag.Number"]
        )

    def saap_validate(
        self, val_evidence: pd.DataFrame, val_msms: pd.DataFrame, saap: pd.DataFrame
    ):
        val_evidence = val_evidence[
            (val_evidence["PEP"] <= self.validation_pep)
            & (val_evidence["PIF"] >= self.pif)
        ].copy()

        if "Reporter intensity corrected 1" in val_evidence.columns:
            norm_reporters = val_evidence.filter(
                regex="Reporter intensity corrected"
            ).apply(lambda x: x / x.median(), axis=0)
            norm_reporters.columns = [
                "Normalised " + col for col in norm_reporters.columns
            ]
            val_evidence = pd.concat([val_evidence, norm_reporters], axis=1)
        else:
            val_evidence.loc[:, "Intensity.Normalised"] = (
                val_evidence["Intensity"] / val_evidence["Intensity"].median()
            )

        filter_list = saap[["Raw file", "Charge", "DP Base Sequence", "SAAP sequence"]]
        filter_list = (
            filter_list.melt(
                id_vars=["Raw file", "Charge"], var_name="Type", value_name="Sequence"
            )
            .drop("Type", axis=1)
            .drop_duplicates()
        )

        val_evidence = val_evidence.merge(
            filter_list, how="inner", on=["Raw file", "Charge", "Sequence"]
        )
        val_evidence = val_evidence.drop(
            ["m/z", "Mass", "Mass error [ppm]", "Retention time", "PIF", "PEP"],
            axis=1,
            errors="ignore",
        )

        msms_ev_merge = val_msms.merge(
            val_evidence, how="inner", on=["Raw file", "MS/MS scan number"]
        )

        saap.loc[:, "b_ions_aas"], saap.loc[:, "y_ions_aas"] = zip(
            *saap.apply(self.frags_containing_aas, axis=1), strict=True
        )
        saap.loc[:, "fragment_evidence"] = self.fragment_evidence(saap, msms_ev_merge)

        # Group and sum fragment evidence across the same SAAP
        # (ignore scan/raw/charge-level variation)
        saap_grouped = (
            saap.groupby(["DP Base Sequence", "SAAP sequence", "aa subs"])[
                "fragment_evidence"
            ]
            .sum()
            .reset_index()
        )
        saap_grouped = saap_grouped[
            saap_grouped["fragment_evidence"] > self.frag_evidence
        ]

        return saap_grouped, msms_ev_merge
