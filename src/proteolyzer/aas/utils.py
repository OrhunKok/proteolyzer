"""Reference tables and small helpers shared by the AAS pipeline stages."""

from pathlib import Path

import numpy as np
import pandas as pd

from proteolyzer import reference
from proteolyzer.core.io import write_frame


def column_mapping(df: pd.DataFrame, cols2keep: list) -> pd.DataFrame:
    """Standardized column selection + renaming based on 'cols2keep'.

    Matching is case-insensitive, and the output column order follows
    ``cols2keep`` rather than the order of the input file.
    """
    col_lookup = {col.lower(): col for col in df.columns}
    rename_map = {
        col_lookup[col.lower()]: col for col in cols2keep if col.lower() in col_lookup
    }

    df = df[list(rename_map.keys())].rename(columns=rename_map)
    df = df[[col for col in cols2keep if col in df.columns]]
    return df


def aa_subs_ref() -> dict:
    """Mass delta of every single amino acid substitution, keyed by origin residue.

    ``{"A": {"A to G": -14.0157, ...}, ...}``
    """
    aa_subs = reference.modifications()
    aa_subs = aa_subs[aa_subs["classification"] == "AA substitution"].copy()
    aa_subs["code_name"] = aa_subs["code_name"].str.replace("2", "->")
    aa_subs["sub_aa"] = (
        aa_subs["code_name"].str.split(" ").str[0].str.split("->").str[1]
    )
    aa_subs["sub_aa"] = aa_subs["sub_aa"].map(reference.three_letter_to_one())
    aa_subs["compact_name"] = aa_subs["one_letter"] + " to " + aa_subs["sub_aa"]
    return (
        aa_subs[["one_letter", "compact_name", "mono_mass"]]
        .groupby("one_letter")
        .apply(
            lambda x: x.set_index("compact_name")["mono_mass"].to_dict(),
            include_groups=False,
        )
        .to_dict()
    )


def gen_mod_dict() -> dict:
    """Known modifications per residue as ``[full_name, position, mono_mass]`` rows.

    Excludes the ``AA substitution`` class, which is what :func:`aa_subs_ref`
    returns. This reference exists to decide whether a mass shift has an
    explanation *other than* mistranslation, so leaving the substitutions in it
    would make every candidate substitution match itself as a modification.
    """
    mod_df = reference.modifications()
    mod_df = mod_df[mod_df["classification"] != "AA substitution"]

    return {
        residue: group[["full_name", "position", "mono_mass"]].values.tolist()
        for residue, group in mod_df.groupby("one_letter")
    }


def ptm_mtp_output(
    dp_df: pd.DataFrame,
    sample_name: str,
    output_dir: Path,
) -> None:
    """Split detected dependent peptides into PTM and MTP sets and pickle both.

    MTP are peptides with a potential AA substitution that cannot be explained
    by a PTM and that have no homologous sequence in any translated frame.
    """
    ptm_df = dp_df[dp_df["PTM"].notna()]
    write_frame(ptm_df, output_dir / "PTM" / f"{sample_name}_PTM")

    frame_cols = [f"{frame}-frame genome substring" for frame in range(1, 7)]
    mtp_df = dp_df[dp_df["PTM"].isna() & dp_df["mistranslated sequence"].notna()]
    mtp_df = mtp_df[~mtp_df[frame_cols].any(axis=1)]
    write_frame(mtp_df, output_dir / "MTP" / f"{sample_name}_MTP")


def calculate_aa_substitution_matrix(
    processed_amino_acids_df: pd.DataFrame,
) -> pd.DataFrame:
    """
    Calculates the pairwise mass difference matrix (Row AA mass - Column AA mass).
    """
    aa_vector: pd.Series = processed_amino_acids_df.set_index("one_letter")["mono_mass"]

    aa_subs_pairwise: pd.DataFrame = pd.DataFrame(
        aa_vector.values[np.newaxis, :] - aa_vector.values[:, np.newaxis],
        index=aa_vector.index.values,
        columns=aa_vector.index.values,
    )

    return aa_subs_pairwise
