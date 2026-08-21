"""Quantification of validated substitutions relative to their base peptide."""

import numpy as np
import pandas as pd

from proteolyzer.core.io import frame_exists, read_frame, write_frame

from .base import Stage


class Quantification(Stage):
    """Computes mistranslated-to-base-peptide intensity ratios per sample."""

    def __init__(self, params, queue=None):
        super().__init__(params, queue)
        # A ratio is only meaningful if both peptides carry enough signal, so
        # the threshold is applied to the mistranslated and the base peptide.
        self.min_quant = float(self.params["Quantification"]["Minimum Quantity"])

    def process_sample(self, sample):
        sample_val_dir = self._locate_sample_dir(sample, suffix="_val")
        if not sample_val_dir:
            self.queue.put(("stderr", f"{sample} validation directory not found"))
            return

        val_evidence_path = (
            self.output_dir / "MTP" / f"{sample}_Val_Evidence_Filtered_Stage_2"
        )
        mtp_path = self.output_dir / "MTP" / f"{sample}_MTP_Filtered_Stage_2"

        if not all(frame_exists(path) for path in (val_evidence_path, mtp_path)):
            self.queue.put(("stderr", f"Missing files for sample {sample}"))
            return

        val_evidence = read_frame(val_evidence_path)
        mtp = read_frame(mtp_path)

        ev_filter_seqs = np.unique(
            mtp[["DP Base Sequence", "mistranslated sequence"]].to_numpy()
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

        mtp = mtp[
            (mtp["DP Base Sequence"].isin(val_evidence.index))
            & (mtp["mistranslated sequence"].isin(val_evidence.index))
        ]
        mtp_df = val_evidence.loc[mtp["mistranslated sequence"]]
        bp_df = val_evidence.loc[mtp["DP Base Sequence"]]

        label_designation = (
            "TMT"
            if "Reporter intensity corrected 1" in val_evidence.columns
            else "Label-Free"
        )
        with np.errstate(divide="ignore"):
            quant = self._raas(
                mtp, mtp_df, bp_df, self.metadata, label_designation=label_designation
            )

        quant = quant[np.isfinite(quant["Ratio"])]
        quant = self._apply_minimum_quantity(quant)

        write_frame(quant, self.output_dir / "MTP" / f"{sample}_MTP_Quant")

        self.queue.put(("stdout", f"Quantification complete for sample: {sample}"))

    def _apply_minimum_quantity(self, quant: pd.DataFrame) -> pd.DataFrame:
        """Drop ratios where either peptide is below ``Minimum Quantity``."""
        if not self.min_quant:
            return quant

        keep = (quant["MTP.Sum"] >= self.min_quant) & (
            quant["BP.Sum"] >= self.min_quant
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

    def _raas(self, mtp, mtp_df, bp_df, sample_df, label_designation):
        """Relative abundance of the substituted peptide vs its base peptide."""
        output_dict = {
            "DP Base Sequence": mtp["DP Base Sequence"],
            "mistranslated sequence": mtp["mistranslated sequence"],
            "aa subs": mtp["aa subs"],
        }

        if label_designation == "Label-Free":
            ratios = np.log2(mtp_df["Intensity"].values / bp_df["Intensity"].values)
        elif label_designation == "TMT":
            ratios = np.log10(mtp_df["Intensity"].values / bp_df["Intensity"].values)

            reporter_regex = r"^(?!.*Normalised).*Reporter intensity corrected.*$"
            mtp_reporters = mtp_df.filter(regex=reporter_regex, axis=1)
            bp_reporters = bp_df.filter(regex=reporter_regex, axis=1)

            norm_regex = "Normalised Reporter intensity corrected"
            mtp_reporters_norm = mtp_df.filter(regex=norm_regex, axis=1)
            bp_reporters_norm = bp_df.filter(regex=norm_regex, axis=1)

            mtp_ratios = mtp_reporters.div(mtp_reporters.sum(axis=1).values, axis=0)
            bp_ratios = bp_reporters.div(bp_reporters.sum(axis=1).values, axis=0)

            mtp_distributed = mtp_ratios.mul(mtp_df["Intensity"], axis=0)
            bp_distributed = bp_ratios.mul(bp_df["Intensity"], axis=0)

            for tmt_plex in sample_df["MQ"].dropna().unique():
                channel = str(int(tmt_plex))
                reporter = f"Reporter intensity corrected {channel}"
                output_dict[f"MTP.Plex.{channel}"] = mtp_distributed[reporter].values
                output_dict[f"BP.Plex.{channel}"] = bp_distributed[reporter].values
                output_dict[f"Ratio.Plex.{channel}"] = np.log10(
                    mtp_distributed[reporter].values / bp_distributed[reporter].values
                )
                output_dict[f"MTP.Plex.{channel}.Norm.Sum"] = mtp_reporters_norm[
                    f"Normalised {reporter}"
                ].values
                output_dict[f"BP.Plex.{channel}.Norm.Sum"] = bp_reporters_norm[
                    f"Normalised {reporter}"
                ].values
        else:
            raise ValueError(f"Unknown label designation: {label_designation!r}")

        output_dict["MTP.Sum"] = mtp_df["Intensity"].values
        output_dict["BP.Sum"] = bp_df["Intensity"].values
        output_dict["Ratio"] = ratios

        return pd.DataFrame.from_dict(output_dict)
