"""Configuration for the AAS pipeline: input files, columns and label schemes."""

from dataclasses import dataclass, field
from typing import ClassVar

from proteolyzer import config as CONFIG


@dataclass(frozen=True)
class DIANN:
    FILES: list[str] = field(default_factory=lambda: ["report", "report-first-pass"])
    #: Columns to load per file. ``None`` loads every column.
    LOAD_COLS: dict[str, list[str] | None] = field(
        default_factory=lambda: {
            "report.stats": None,
            "report": None,
            "report-first-pass": None,
        }
    )
    FILE_EXTENSIONS: list[str] = field(default_factory=lambda: [".parquet", ".tsv"])
    COLS_RENAME_MAPPING: dict[str, str] = field(default_factory=dict)
    #: Files each stage needs to find in a search output directory.
    FILES_NEEDED: dict[str, list[str]] = field(
        default_factory=lambda: {
            "Detection": ["report"],
            "Validation": ["report"],
        }
    )


@dataclass(frozen=True)
class MaxQuant:
    FILES: list[str] = field(
        default_factory=lambda: ["allPeptides", "evidence", "msms", "peptides"]
    )
    #: Columns to load per file. ``None`` loads every column.
    LOAD_COLS: dict[str, list[str] | None] = field(
        default_factory=lambda: {
            "allPeptides": [
                "Raw file",
                "Charge",
                "m/z",
                "Mass",
                "Retention time",
                "Intensity",
                "DP Mass Difference",
                "DP PEP",
                "DP Decoy",
                "DP Proteins",
                "Reverse",
                "DP Base Sequence",
                "DP Probabilities",
                "DP Positional Probability",
                "DP Base Scan Number",
                "DP Mod Scan Number",
                "MSMS Scan Numbers",
            ],
            "evidence": [
                "Raw file",
                "Charge",
                "m/z",
                "Mass",
                "Retention time",
                "Reverse",
                "Potential contaminant",
                "Sequence",
                "PIF",
                "PEP",
                "Mass error [ppm]",
                "MS/MS scan number",
                "Intensity",
            ],
            "msms": ["Raw file", "Scan number", "Matches", "Reverse"],
            "peptides": [
                "Sequence",
                "Start position",
                "Amino acid after",
                "Amino acid before",
                "Reverse",
                "Potential contaminant",
            ],
        }
    )
    FILE_EXTENSIONS: list[str] = field(default_factory=lambda: [".txt"])
    COLS_RENAME_MAPPING: dict[str, str] = field(
        default_factory=lambda: {
            "Experiment": "Run",
            "Modified sequence": "Precursor.Id",
            "Sequence": "Stripped.Sequence",
            "Charge": "Precursor.Charge",
            "Gene names": "Genes",
            "Leading razor protein": "Leading.Razor.Protein",
            "Missed cleavages": "Missed.Cleavages",
            "Retention time": "RT",
            "Retention length": "RT.Width",
        }
    )
    #: Files each stage needs to find in a search output directory.
    FILES_NEEDED: dict[str, list[str]] = field(
        default_factory=lambda: {
            "Detection": ["evidence", "allPeptides", "peptides"],
            "Validation": ["evidence", "msms"],
        }
    )


@dataclass(frozen=True)
class TMT:
    """TMT reference tables.

    Only ``MQ_TMT_MAP`` is consumed at the moment (by ``Stage.metadata``);
    the other two are kept as reference data for TMT-specific workflows.
    """

    TMT_PLEX_MAP: ClassVar[dict[str, str]] = {
        "126": "1", "127N": "2", "127C": "3", "128N": "4", "128C": "5",
        "129N": "6", "129C": "7", "130N": "8", "130C": "9", "131N": "10",
        "131C": "11", "132N": "12", "132C": "13", "133N": "14", "133C": "15",
        "134N": "16", "134C": "17", "135N": "18",
    }  # fmt: skip

    MQ_TMT_MAP: ClassVar[dict[str, int]] = {
        "126": 1, "127N": 2, "127C": 3, "128N": 4, "128C": 5,
        "129N": 6, "129C": 7, "130N": 8, "130C": 9, "131": 10,
    }  # fmt: skip

    MASS_SHIFT: ClassVar[dict[int, float]] = {
        0: 224.152478,
        2: 225.155833,
        6: 229.162932,
        10: 229.162932,
        11: 229.169252,
    }


@dataclass(frozen=True)
class Config:
    DIANN: DIANN = field(default_factory=DIANN)
    MaxQuant: MaxQuant = field(default_factory=MaxQuant)
    TMT: TMT = field(default_factory=TMT)
    Protease: CONFIG.Protease = field(default_factory=CONFIG.Protease)
    AminoAcids: CONFIG.AminoAcids = field(default_factory=CONFIG.AminoAcids)
