"""Input-format configuration for the core loading/processing pipeline.

Each search-engine block describes the files proteolyzer recognizes, the
columns worth loading from them, how those columns map onto the canonical
proteolyzer names, and which of them must stay numeric.
"""

from dataclasses import dataclass, field

#: Columns loaded from the DIA-NN precursor reports. The main report and the
#: first-pass report share a schema, so they share this set.
_DIANN_REPORT_COLS: frozenset[str] = frozenset(
    {
        "Run",
        "Precursor.Id",
        "Stripped.Sequence",
        "Precursor.Charge",
        "Genes",
        "Protein.Group",
        "Proteotypic",
        "RT",
        "RT.Start",
        "RT.Stop",
        "FWHM",
        "IM",
        "PEP",
        "Ms1.Area",
        "Ms1.Apex.Area",
        "Ms1.Normalised",
        "Precursor.Quantity",
        "Precursor.Normalised",
    }
)


@dataclass(frozen=True)
class DIANN:
    FILES: list[str] = field(
        default_factory=lambda: [
            "report",
            "report-first-pass",
            "report-first-pass.site_report",
            "report.site_report",
            "report.stats",
            "report.stats-first-pass",
            "report.log",
            "xic",
        ]
    )
    #: Columns to load per file. ``None`` loads every column.
    LOAD_COLS: dict[str, set | None] = field(
        default_factory=lambda: {
            "report.stats": None,
            "report": set(_DIANN_REPORT_COLS),
            "report-first-pass": set(_DIANN_REPORT_COLS),
            "report-first-pass.site_report": None,
            "report.site_report": None,
            "report.log": None,
            "xic": None,
        }
    )
    FILE_EXTENSIONS: list[str] = field(default_factory=lambda: [".parquet", ".tsv"])
    COLS_RENAME_MAPPING: dict[str, str] = field(default_factory=dict)
    EXCLUDE_CAT_CONVERSION: set[str] = field(
        default_factory=lambda: {
            "Precursor.Quantity",
            "Ms1.Apex.Area",
            "Ms1.Normalised",
            "Ms1.Area",
            "Precursor.Normalised",
        }
    )


@dataclass(frozen=True)
class MaxQuant:
    FILES: list[str] = field(
        default_factory=lambda: [
            "allPeptides",
            "evidence",
            "dependentPeptides",
            "matchedFeatures",
            "modificationSpecificPeptides",
            "ms3Scans",
            "msms",
            "msmsScans",
            "mzRange",
            "Oxidation (M)Sites",
            "parameters",
            "peptides",
            "proteinGroups",
            "summary",
        ]
    )
    #: Columns to load per file. ``None`` loads every column.
    LOAD_COLS: dict[str, set | None] = field(
        default_factory=lambda: {
            "allPeptides": None,
            "evidence": {
                "Experiment",
                "Modified sequence",
                "Sequence",
                "Charge",
                "Gene names",
                "Proteins",
                "Leading razor protein",
                "Missed cleavages",
                "Retention time",
                "Retention length",
                "PEP",
                "Intensity",
            },
            "dependentPeptides": None,
            "matchedFeatures": None,
            "modificationSpecificPeptides": None,
            "ms3Scans": None,
            "msms": None,
            "msScans": None,
            "msmsScans": None,
            "mzRange": None,
            "Oxidation (M)Sites": None,
            "parameters": None,
            "peptides": None,
            "proteinGroups": None,
            "summary": None,
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
    EXCLUDE_CAT_CONVERSION: set[str] = field(default_factory=lambda: {"Intensity"})


@dataclass(frozen=True)
class Config:
    COL_MEDIAN_THRESHOLD: int = 100
    CARDINALITY_THRESHOLD: float = 0.1
    DIANN: DIANN = field(default_factory=DIANN)
    MaxQuant: MaxQuant = field(default_factory=MaxQuant)
