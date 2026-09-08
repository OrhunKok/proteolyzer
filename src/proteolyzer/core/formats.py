"""Input-format configuration for the core loading/processing pipeline.

Each block describes one search engine's output: the files it goes by, the
columns that identify it, how its columns map onto the canonical proteolyzer
names, which of them are written as text and are not text, and any canonical
column it does not write that can be built out of ones it does. All of that is a
fact about the format, true for everyone who reads it.

How a file is matched to a block, and why the columns rather than the name, is in
``docs/notes/recognising-a-format.md``.

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
class Narrow:
    """A table recognized by the shape of its whole schema, not a few names in it.

    For a table too small and too plainly named for an overlap of names to say
    anything. Both conditions have to hold, and why both is in
    ``docs/notes/recognising-a-format.md``.
    """

    #: Every one of these has to be present, not two of them.
    columns: frozenset[str]
    #: The most columns the table is expected to have. Generous: it guards
    #: against something an order of magnitude wider.
    width: int


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
    #: DIA-NN's own names are the canonical schema, so these are the schema's
    #: too. The report only; `xic` is signed by its shape instead, below.
    COLUMN_SIGNATURE: frozenset[str] = field(
        default_factory=lambda: frozenset(
            {
                "Precursor.Id",
                "Precursor.Lib.Index",
                "Precursor.Normalised",
                "Precursor.Quantity",
                "Modified.Sequence",
                "Stripped.Sequence",
                "Proteotypic",
                "Protein.Group",
                "Ms1.Area",
                "Ms1.Normalised",
                "PG.MaxLFQ",
                "PG.Q.Value",
                "Global.Q.Value",
                "Lib.Q.Value",
                "Run.Index",
            }
        )
    )
    #: The XIC export. Its four column names each belong to nobody in
    #: particular, so its shape is what identifies it; the width ceiling is
    #: loose on purpose. See ``docs/notes/recognising-a-format.md``.
    NARROW_SIGNATURES: tuple[Narrow, ...] = field(
        default_factory=lambda: (
            Narrow(columns=frozenset({"pr", "feature", "rt", "value"}), width=8),
        )
    )
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
    #: Fifteen tables that share few columns, so this reaches any of them
    #: rather than describing one: `Raw file` carries most, the scan tables
    #: their instrument readings, `proteinGroups` its protein ones.
    #:
    #: `Charge` and `Intensity` are absent deliberately -- FragPipe writes both
    #: under those names. See ``docs/notes/recognising-a-format.md``.
    COLUMN_SIGNATURE: frozenset[str] = field(
        default_factory=lambda: frozenset(
            {
                "Raw file",
                "Modified sequence",
                "Leading razor protein",
                "MS/MS count",
                "MS scan number",
                "Ion injection time",
                "Total ion current",
                "Precursor apex offset time",
                "Retention length",
                "Retention length (FWHM)",
                "Mass deficit",
                "PIF",
                "Protein IDs",
                "Majority protein IDs",
                "Potential contaminant",
            }
        )
    )
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
    #: All three identification tables share one set of columns. `rt` and `mz`
    #: are left out, being DIA-NN's XIC export's as well.
    COLUMN_SIGNATURE: frozenset[str] = field(
        default_factory=lambda: frozenset(
            {
                "file_name",
                "stripped_seq",
                "untag_seq",
                "untag_prec",
                "silac_channel",
                "channels_matched",
                "plex_Area",
                "scribe_scores",
                "frac_dia_int",
                "BestChannel_Qvalue",
                "window_mz",
                "iso_cor",
                "pep_len",
                "is_decoy",
            }
        )
    )
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
    #: The search scores for `psm`, the counting columns for `protein`, so
    #: either table is reached. `Charge`, `Intensity`, `Peptide`, `Protein` and
    #: `Gene` are absent for being too plain to be anybody's in particular.
    COLUMN_SIGNATURE: frozenset[str] = field(
        default_factory=lambda: frozenset(
            {
                "Spectrum",
                "Spectrum File",
                "Modified Peptide",
                "Assigned Modifications",
                "Hyperscore",
                "Nextscore",
                "PeptideProphet Probability",
                "Calibrated Observed M/Z",
                "Number of Enzymatic Termini",
                "Razor Spectral Count",
                "Total Spectral Count",
                "Razor Intensity",
                "Protein Probability",
                "Top Peptide Probability",
                "Unique Peptides",
            }
        )
    )
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


#: The report's columns onto the canonical schema, spelled the way the
#: tab-separated export spells them. The parquet export writes the same names
#: with ``_`` for ``.``, and :class:`Spectronaut` carries both -- derived from
#: this rather than written twice. See ``docs/notes/spectronaut.md``.
_SPECTRONAUT_COLUMNS: dict[str, str] = {
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
    # PG.Quantity is absent because it is already the canonical name: the
    # schema is DIA-NN's own, and the two agree on this one.
    "PG.Qvalue": "PG.Q.Value",
    "PEP.NrOfMissedCleavages": "Missed.Cleavages",
    # True and False, where DIA-NN writes 1 and 0. The name is the same on both
    # sides of the rename and the dtype is not, so a caller comparing one
    # against 0 has to say `== False` instead.
    "PEP.IsProteotypic": "Proteotypic",
    "EG.IsDecoy": "Decoy",
}


def _either_spelling(names: set[str]) -> frozenset[str]:
    """`names` as the text export spells them and as the parquet export does."""
    return frozenset(names) | {name.replace(".", "_") for name in names}


def _both_spellings(mapping: dict[str, str]) -> dict[str, str]:
    """`mapping` keyed by the text export's names and by the parquet's.

    A name the file does not carry does nothing, so holding both costs a dict
    twice the size and buys not having to know which serialization this is.
    """
    return mapping | {
        name.replace(".", "_"): canonical for name, canonical in mapping.items()
    }


@dataclass(frozen=True)
class Spectronaut:
    """A Spectronaut report: long format, one row a precursor a run.

    Columns are prefixed by the level they belong to -- ``E.`` the experiment,
    ``R.`` a run, ``PG.`` a protein group, ``PEP.`` a peptide, ``EG.`` an
    elution group, ``FG.`` a fragment group, which is a precursor.

    Read as parquet or as tab-separated text; the extension decides, and a
    caller does nothing to pick. The two are one report in two serializations
    and come back as the same frame, which a test asserts -- but they do not
    spell their column names the same way, hence both spellings below.

    A report is configurable column by column, so a subset naming a column the
    analysis did not write must not fail the read. The intersection
    ``cols_to_load`` already takes is what makes that so.

    ``docs/notes/spectronaut.md`` has what was measured off real exports,
    including the quantity scale, the isolation-window column, and the two
    things this block had wrong before v0.9.0.
    """

    #: A name the export sometimes goes by, not a name it must have. Kept
    #: because a report exported as one is then recognized without opening it.
    FILES: list[str] = field(default_factory=lambda: ["Report"])
    #: A shape an export often has and never has to: Spectronaut has no
    #: default output name, so this is a shortcut and COLUMN_SIGNATURE is the
    #: identification. Matched against the stem in full, and case-sensitively.
    FILE_PATTERNS: list[str] = field(default_factory=lambda: [r".*_Report"])
    #: Both. DIA-NN claims ``.parquet`` too, which the case-sensitivity above
    #: is what holds apart.
    FILE_EXTENSIONS: list[str] = field(default_factory=lambda: [".parquet", ".tsv"])
    #: What identifies the report when its name cannot, which is the usual
    #: case. Any two settle it: level-prefixed the way no other engine here
    #: prefixes anything, in either spelling.
    COLUMN_SIGNATURE: frozenset[str] = field(
        default_factory=lambda: frozenset(
            _both_spellings(
                {
                    name: name
                    for name in (
                        "R.FileName",
                        "EG.ModifiedSequence",
                        "PEP.StrippedSequence",
                        "FG.Charge",
                        "FG.Quantity",
                        "PG.ProteinGroups",
                    )
                }
            )
        )
    )
    COLS_RENAME_MAPPING: dict[str, str] = field(
        default_factory=lambda: _both_spellings(_SPECTRONAUT_COLUMNS)
    )
    #: No EG.PrecursorId in every export, so the identifier the rest of the
    #: package keys on is built from the two columns a fragment group always
    #: has. Under the canonical names, being built after the rename.
    BUILT_COLS: dict[str, tuple[str, ...]] = field(
        default_factory=lambda: {
            "Precursor.Id": ("Modified.Sequence", "Precursor.Charge")
        }
    )
    #: Columns the export writes as text that are not text, named as the file
    #: names them so a caller keeping the file's own names gets them too. A
    #: list rather than a rule, because a rule turns the XIC database key into
    #: an integer; see ``docs/notes/numbers-and-gaps.md``.
    NUMERIC_COLS: frozenset[str] = field(
        default_factory=lambda: _either_spelling(
            {
                "EG.Qvalue",
                "EG.MinChannelQvalue",
                "EG.MaxChannelQvalue",
                "PEP.NrOfMissedCleavages",
            }
        )
    )
    #: As above. The tab-separated export writes these as real booleans, so
    #: this is also what keeps the two serializations agreeing.
    BOOLEAN_COLS: frozenset[str] = field(
        default_factory=lambda: _either_spelling({"PEP.IsProteotypic"})
    )
    #: Nothing: the columns worth keeping out are numbers, which are never
    #: converted whatever this says.
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
