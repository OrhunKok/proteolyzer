"""Quantification of validated substitutions relative to their base peptide."""

import numpy as np
import pandas as pd

from proteolyzer.core.io import frame_exists, read_frame, write_frame

from .base import Stage


class Quantification(Stage):
    """Computes SAAP-to-base-peptide intensity ratios per sample."""

    def __init__(self, params, queue=None):
        super().__init__(params, queue)
        # A ratio is only meaningful if both peptides carry enough signal, so
        # the threshold is applied to the SAAP and to the BASE peptide.
        self.min_quant = float(self.params["Quantification"]["Minimum Quantity"])

    def process_sample(self, sample):
        sample_val_dir = self._locate_sample_dir(sample, suffix="_val")
        if not sample_val_dir:
            self.queue.put(("stderr", f"{sample} validation directory not found"))
            return

        val_evidence_path = (
            self.output_dir / "SAAP" / f"{sample}_Val_Evidence_Filtered_Stage_2"
        )
        saap_path = self.output_dir / "SAAP" / f"{sample}_SAAP_Filtered_Stage_2"

        if not all(frame_exists(path) for path in (val_evidence_path, saap_path)):
            self.queue.put(("stderr", f"Missing files for sample {sample}"))
            return

        val_evidence = read_frame(val_evidence_path)
        saap = read_frame(saap_path)

        ev_filter_seqs = np.unique(
            saap[["DP Base Sequence", "SAAP sequence"]].to_numpy()
        )
        val_evidence = (
            val_evidence[val_evidence["Sequence"].isin(ev_filter_seqs)]
            .drop(["Frag.Type", "Frag.Number"], axis=1)
            .drop_duplicates()
        )

        val_evidence = (
            val_evidence.drop(["Raw file", "Charge", "MS/MS scan number"], axis=1)
            .groupby("Sequence")
            .sum()
        )

        saap = saap[
            (saap["DP Base Sequence"].isin(val_evidence.index))
            & (saap["SAAP sequence"].isin(val_evidence.index))
        ]
        saap_df = val_evidence.loc[saap["SAAP sequence"]]
        base_df = val_evidence.loc[saap["DP Base Sequence"]]

        label_designation = (
            "TMT"
            if "Reporter intensity corrected 1" in val_evidence.columns
            else "Label-Free"
        )
        with np.errstate(divide="ignore"):
            quant = self._raas(
                saap,
                saap_df,
                base_df,
                self.metadata,
                label_designation=label_designation,
            )

        quant = quant[np.isfinite(quant["Ratio"])]
        quant = self._apply_minimum_quantity(quant)

        write_frame(quant, self.output_dir / "SAAP" / f"{sample}_SAAP_Quant")

        self.queue.put(("stdout", f"Quantification complete for sample: {sample}"))

    def _apply_minimum_quantity(self, quant: pd.DataFrame) -> pd.DataFrame:
        """Drop ratios where either peptide is below ``Minimum Quantity``."""
        if not self.min_quant:
            return quant

        keep = (quant["SAAP.Sum"] >= self.min_quant) & (
            quant["BASE.Sum"] >= self.min_quant
        )
        dropped = int((~keep).sum())
        if dropped:
            self.queue.put(
                (
                    "stdout",
                    f"Dropped {dropped} ratios below the minimum quantity of "
                    f"{self.min_quant:g}.",
                )
            )
        return quant[keep]

    def _raas(self, saap, saap_df, base_df, sample_df, label_designation):
        """Relative abundance of the substituted peptide vs its base peptide."""
        output_dict = {
            "DP Base Sequence": saap["DP Base Sequence"],
            "SAAP sequence": saap["SAAP sequence"],
            "aa subs": saap["aa subs"],
        }

        if label_designation == "Label-Free":
            ratios = np.log2(saap_df["Intensity"].values / base_df["Intensity"].values)
        elif label_designation == "TMT":
            ratios = np.log10(saap_df["Intensity"].values / base_df["Intensity"].values)

            reporter_regex = r"^(?!.*Normalised).*Reporter intensity corrected.*$"
            saap_reporters = saap_df.filter(regex=reporter_regex, axis=1)
            base_reporters = base_df.filter(regex=reporter_regex, axis=1)

            norm_regex = "Normalised Reporter intensity corrected"
            saap_reporters_norm = saap_df.filter(regex=norm_regex, axis=1)
            base_reporters_norm = base_df.filter(regex=norm_regex, axis=1)

            saap_ratios = saap_reporters.div(saap_reporters.sum(axis=1).values, axis=0)
            base_ratios = base_reporters.div(base_reporters.sum(axis=1).values, axis=0)

            saap_distributed = saap_ratios.mul(saap_df["Intensity"], axis=0)
            base_distributed = base_ratios.mul(base_df["Intensity"], axis=0)

            for tmt_plex in sample_df["MQ"].dropna().unique():
                channel = str(int(tmt_plex))
                reporter = f"Reporter intensity corrected {channel}"
                output_dict[f"SAAP.Plex.{channel}"] = saap_distributed[reporter].values
                output_dict[f"BASE.Plex.{channel}"] = base_distributed[reporter].values
                output_dict[f"Ratio.Plex.{channel}"] = np.log10(
                    saap_distributed[reporter].values
                    / base_distributed[reporter].values
                )
                output_dict[f"SAAP.Plex.{channel}.Norm.Sum"] = saap_reporters_norm[
                    f"Normalised {reporter}"
                ].values
                output_dict[f"BASE.Plex.{channel}.Norm.Sum"] = base_reporters_norm[
                    f"Normalised {reporter}"
                ].values
        else:
            raise ValueError(f"Unknown label designation: {label_designation!r}")

        output_dict["SAAP.Sum"] = saap_df["Intensity"].values
        output_dict["BASE.Sum"] = base_df["Intensity"].values
        output_dict["Ratio"] = ratios

        return pd.DataFrame.from_dict(output_dict)
