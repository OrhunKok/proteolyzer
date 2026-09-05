"""Input-format configuration for the core loading/processing pipeline.

Each search-engine block describes the files proteolyzer recognizes -- by name,
or by pattern where the engine stamps the name with the moment it wrote it --
how their columns map onto the canonical proteolyzer names, which of them must
stay numeric, and any canonical column the format does not write that can be
built out of ones it does. All of that is a fact about the format, true for
everyone who reads it.

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
    """A Spectronaut report: long format, one row a precursor a run.

    Columns are prefixed by the level they belong to -- ``E.`` the experiment,
    ``R.`` a run, ``PG.`` a protein group, ``PEP.`` a peptide, ``EG.`` an
    elution group, ``FG.`` a fragment group, which is a precursor.

    A report is configurable column by column, so what one lab's export holds is
    not what another's does: 78 columns in the one this was written from. That
    makes the intersection :meth:`~proteolyzer.core.loader.DataLoader._cols_to_load`
    already takes load-bearing rather than convenient -- naming a column the
    analysis did not write must not fail the read.

    Parquet and tab-separated text, because Spectronaut writes either and
    parquet is what it writes by default. The two are one report in two
    serializations rather than two formats: the same names, the same rename
    mapping, the same built identifier, and a test asserting a report read both
    ways comes back the same. Which is read is decided by the extension, so a
    caller does nothing to pick.

    Two things measured off that export -- 13 runs, 173,443 rows, 174 MB, tab
    separated -- are worth knowing before reading one:

    ``FG.Quantity`` spans 2.54 to 400,000 in the one column, so it is not on the
    scale a DIA-NN area is. ``round_large_floats`` would throw away a fifth of a
    precursor quantified at 2.54, and must stay off for this format, as it is by
    default for every format.

    ``FG.PrecWindowNumber`` says which of the method's isolation windows took the
    precursor, which no other report read here states -- grouping by it recovers
    the window scheme with no design file to hand. It is an integer, and a number
    is never made categorical, so nothing has to keep it out of that.
    """

    #: The table, under the name it goes by when Spectronaut has not stamped it.
    FILES: list[str] = field(default_factory=lambda: ["Report"])
    #: Spectronaut writes ``<date>_<time>_<analysis>_Report.tsv``, so no fixed
    #: name can match a real export. Matched against the stem in full rather
    #: than from its start, so the ``..._Report.setup`` written beside it is not
    #: taken for the report; and case-sensitively, because DIA-NN's
    #: ``report.tsv`` differs from a bare ``Report.tsv`` by the one letter, and
    #: a file two blocks claim is an error rather than a guess.
    FILE_PATTERNS: list[str] = field(default_factory=lambda: [r".*_Report"])
    #: Both, and the pattern is over the stem so it does not care which. The
    #: case-sensitivity above carries more weight now than it did: DIA-NN
    #: claims ``.parquet`` too, and its ``report.parquet`` is one capital
    #: letter from a bare ``Report.parquet``.
    FILE_EXTENSIONS: list[str] = field(default_factory=lambda: [".parquet", ".tsv"])
    COLS_RENAME_MAPPING: dict[str, str] = field(
        default_factory=lambda: {
            "R.FileName": "Run",
            "EG.ModifiedSequence": "Modified.Sequence",
            "PEP.StrippedSequence": "Stripped.Sequence",
            "FG.Charge": "Precursor.Charge",
            "FG.PrecMz": "Precursor.Mz",
            "FG.Quantity": "Precursor.Quantity",
            "EG.ApexRT": "RT",
            "EG.RTPredicted": "Predicted.RT",
            "EG.iRTEmpirical": "iRT",
            "EG.Qvalue": "Q.Value",
            "EG.PEP": "PEP",
            "PG.ProteinGroups": "Protein.Group",
            "PG.ProteinAccessions": "Protein.Ids",
            # PG.Quantity is absent because it is already the canonical name:
            # the schema is DIA-NN's own, and the two agree on this one.
            "PG.Qvalue": "PG.Q.Value",
            "PEP.NrOfMissedCleavages": "Missed.Cleavages",
            # True and False, where DIA-NN writes 1 and 0. The name is the same
            # on both sides of the rename and the dtype is not, so a caller
            # comparing one against 0 has to say `== False` instead.
            "PEP.IsProteotypic": "Proteotypic",
            "EG.IsDecoy": "Decoy",
        }
    )
    #: There is no EG.PrecursorId in every export -- there was none in the one
    #: this was written from -- so the identifier the rest of the package keys
    #: on is built from the two columns a fragment group always has. Stated
    #: under the canonical names, since it is built after the rename.
    BUILT_COLS: dict[str, tuple[str, ...]] = field(
        default_factory=lambda: {
            "Precursor.Id": ("Modified.Sequence", "Precursor.Charge")
        }
    )
    #: Nothing. The columns worth keeping out are the quantitative ones, and
    #: they are numbers, which are never converted whatever this says.
    EXCLUDE_CAT_CONVERSION: set[str] = field(default_factory=set)


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
