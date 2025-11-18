from dataclasses import dataclass, field
from typing import List, Dict
from .. import config as CONFIG

@dataclass(frozen=True) 
class DIANN:
    FILES: List[str] = field(default_factory=lambda: [
        "report", "report-first-pass", "report.stats", "report.stats-first-pass", "report.log"
    ])
    LOAD_COLS: Dict[str, set] = field(default_factory=lambda: {
        "report.stats": None,  # Load all columns
        "report": {
            "Run", "Precursor.Id", "Stripped.Sequence", "Precursor.Charge", "Genes",
            "Protein.Group", "Proteotypic", "RT", "RT.Start", "RT.Stop", "FWHM", "IM",
            "PEP", "Ms1.Area", "Ms1.Apex.Area", "Ms1.Normalised", "Precursor.Quantity",
            "Precursor.Normalised"
        },
        "report-first-pass": {
            "Run", "Precursor.Id", "Stripped.Sequence", "Precursor.Charge", "Genes",
            "Protein.Group", "Proteotypic", "RT", "RT.Start", "RT.Stop", "FWHM", "IM",
            "PEP", "Ms1.Area", "Ms1.Apex.Area", "Ms1.Normalised", "Precursor.Quantity",
            "Precursor.Normalised"
        },
        "xic": None,  # Load all columns
    })
    FILE_EXTENSIONS: List[str] = field(default_factory=lambda: [".parquet", ".tsv"])
    COLS_RENAME_MAPPING: Dict[str, str] = field(default_factory=lambda: {})
    EXCLUDE_CAT_CONVERSION: set = field(default_factory=lambda: {
        'Precursor.Quantity', 'Ms1.Apex.Area', 'Ms1.Normalised', 'Ms1.Area', 'Precursor.Normalised'
    })

@dataclass(frozen=True)
class MaxQuant:
    FILES: List[str] = field(default_factory=lambda: [
        "allPeptides", "evidence", "matchedFeatures", "modificationSpecificPeptides",
        "ms3Scans", "msms", "msmsScans", "mzRange", "Oxidation (M)Sites", "parameters",
        "peptides", "proteinGroups", "summary"
    ])
    LOAD_COLS: Dict[str, set] = field(default_factory=lambda: {
        "allPeptides": None,  # Load all columns
        "evidence": {
            "Experiment", "Modified sequence", "Sequence", "Charge", "Gene names",
            "Proteins", "Leading razor protein", "Missed cleavages", "Retention time",
            "Retention length", "PEP", "Intensity"
        },
        "msScans": None,  # Load all columns
        "msmsScans": None,  # Load all columns
        "parameters": None,  # Load all columns
        "summary": None,  # Load all columns
    })
    FILE_EXTENSIONS: List[str] = field(default_factory=lambda: [".txt"])
    COLS_RENAME_MAPPING: Dict[str, str] = field(default_factory=lambda: {
        "Experiment": "Run",
        "Modified sequence": "Precursor.Id",
        "Sequence": "Stripped.Sequence",
        "Charge": "Precursor.Charge",
        "Gene names": "Genes",
        "Protein.Group": "Protein.Group",
        "Leading razor protein": "Leading.Razor.Protein",
        "Missed cleavages": "Missed.Cleavages",
        "Retention time": "RT",
        "Retention length": "RT.Width",
    })
    EXCLUDE_CAT_CONVERSION: set = field(default_factory=lambda: {'Intensity'})


@dataclass(frozen=True)
class Config:
    COL_MEDIAN_THRESHOLD: int = 100
    CARDINALITY_THRESHOLD: float = 0.1
    DIANN: DIANN = field(default_factory=DIANN)
    MaxQuant: MaxQuant = field(default_factory=MaxQuant)
    Protease: CONFIG.Protease = field(default_factory=CONFIG.Protease)