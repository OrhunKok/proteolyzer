"""Input-format configuration for the core loading/processing pipeline.

Each search-engine block describes the files proteolyzer recognizes, how their
columns map onto the canonical proteolyzer names, and which of them must stay
numeric. All of that is a fact about the format, true for everyone who reads it.

**Which columns to keep is not here, deliberately.** It is a fact about the
project doing the reading rather than about the file: a dashboard plots the m/z
and the injection time a pipeline never looks at, and a pipeline wants the
quantities a dashboard has no panel for. Carrying one list upstream made every
consumer either take a subset built for somebody else or override it -- which is
what happened -- so the list belongs to the caller, and ``Data.cols_to_load``
is how it is stated. Reading everything is the default, because a reader that
silently drops a column is worse than a wide frame. See DECISIONS.md.
"""

from dataclasses import dataclass, field


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
            "msScans",
            "msmsScans",
            "mzRange",
            "Oxidation (M)Sites",
            "parameters",
            "peptides",
            "proteinGroups",
            "summary",
        ]
    )
    FILE_EXTENSIONS: list[str] = field(default_factory=lambda: [".txt"])
    COLS_RENAME_MAPPING: dict[str, str] = field(
        default_factory=lambda: {
            "Experiment": "Run",
            #: evidence names the run Experiment; every other table names it
            #: Raw file.
            "Raw file": "Run",
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
class JMod:
    FILES: list[str] = field(
        default_factory=lambda: ["filtered_IDs", "all_IDs", "all_IDs_filtered"]
    )
    FILE_EXTENSIONS: list[str] = field(default_factory=lambda: [".csv", ".parquet"])
    COLS_RENAME_MAPPING: dict[str, str] = field(
        default_factory=lambda: {
            "file_name": "Run",
            "seq": "Precursor.Id",
            "stripped_seq": "Stripped.Sequence",
            "z": "Precursor.Charge",
            "protein": "Protein.Group",
            "silac_channel": "Channel",
            "rt": "RT",
            "mz": "Precursor.Mz",
            "plex_Area": "Precursor.Quantity",
            "MS1_Area": "Ms1.Area",
            "pep_len": "Peptide.Length",
        }
    )
    EXCLUDE_CAT_CONVERSION: set[str] = field(default_factory=set)


@dataclass(frozen=True)
class FragPipe:
    FILES: list[str] = field(
        default_factory=lambda: ["psm", "peptide", "ion", "protein"]
    )
    FILE_EXTENSIONS: list[str] = field(default_factory=lambda: [".tsv"])
    COLS_RENAME_MAPPING: dict[str, str] = field(
        default_factory=lambda: {
            "Spectrum File": "Run",
            "Modified Peptide": "Precursor.Id",
            "Peptide": "Stripped.Sequence",
            "Charge": "Precursor.Charge",
            "Retention": "RT",
            "Observed M/Z": "Precursor.Mz",
            "Intensity": "Precursor.Quantity",
            "Protein": "Protein.Group",
            "Gene": "Genes",
            "Peptide Length": "Peptide.Length",
            "Number of Missed Cleavages": "Missed.Cleavages",
            "Ion Mobility": "IM",
        }
    )
    #: One row per spectrum, so it is unique per row and categorising it costs
    #: memory rather than saving it.
    EXCLUDE_CAT_CONVERSION: set[str] = field(default_factory=lambda: {"Spectrum"})


@dataclass(frozen=True)
class Spectronaut:
    #: Spectronaut names a report ``<date>_<time>_<analysis>_Report.tsv``, so
    #: no exact name can match every export the way FILES matches the other
    #: engines'. FILES stays empty and the file's stem is matched against
    #: FILE_PATTERN instead.
    FILE_PATTERN: str = r".*_?Report$"
    FILES: list[str] = field(default_factory=list)
    FILE_EXTENSIONS: list[str] = field(default_factory=lambda: [".tsv"])
    COLS_RENAME_MAPPING: dict[str, str] = field(
        default_factory=lambda: {
            "R.FileName": "Run",
            "EG.ModifiedSequence": "Modified.Sequence",
            "PEP.StrippedSequence": "Stripped.Sequence",
            "FG.Charge": "Precursor.Charge",
            "FG.PrecMz": "Precursor.Mz",
            "EG.ApexRT": "RT",
            "EG.RTPredicted": "Predicted.RT",
            "EG.iRTEmpirical": "iRT",
            "FG.Quantity": "Precursor.Quantity",
            "PG.ProteinGroups": "Protein.Group",
            "PG.ProteinAccessions": "Protein.Ids",
            "EG.Qvalue": "Q.Value",
            "EG.PEP": "PEP",
            "PG.Qvalue": "PG.Q.Value",
            "PEP.IsProteotypic": "Proteotypic",
            "EG.IsDecoy": "Decoy",
            #: Precursor.Id is not here: this export carries no
            #: EG.PrecursorId, and building one is EG.ModifiedSequence +
            #: FG.Charge -- two columns, not a rename.
        }
    )
    #: Quantities in this format measured 2.54 to 400,000 in one column, and
    #: FG.PrecWindowNumber is the isolation window a precursor was taken
    #: from -- both lose their meaning made categorical.
    EXCLUDE_CAT_CONVERSION: set[str] = field(
        default_factory=lambda: {
            "Precursor.Quantity",
            "PG.Quantity",
            "FG.PrecWindowNumber",
        }
    )


@dataclass(frozen=True)
class Config:
    COL_MEDIAN_THRESHOLD: int = 100
    #: Fraction of a column's memory that turning it categorical has to save
    #: for the conversion to be worth making. Set from measurement: on a real
    #: report the columns that benefit save 49% or more and the ones that do
    #: not save under 2% (or cost memory), so anything in between separates
    #: them.
    MIN_CATEGORICAL_SAVING: float = 0.2
    DIANN: DIANN = field(default_factory=DIANN)
    MaxQuant: MaxQuant = field(default_factory=MaxQuant)
    JMod: JMod = field(default_factory=JMod)
    FragPipe: FragPipe = field(default_factory=FragPipe)
    Spectronaut: Spectronaut = field(default_factory=Spectronaut)
