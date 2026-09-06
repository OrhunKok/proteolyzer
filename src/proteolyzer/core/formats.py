"""Input-format configuration for the core loading/processing pipeline.

Each search-engine block describes the files proteolyzer recognizes -- by name
where the engine names its own output, and by the columns inside where the name
is a convention rather than a fact -- how their columns map onto the canonical
proteolyzer names, which of them are written as text and are not text, and any
canonical column the format does not write that can be built out of ones it
does. All of that is a fact about the format, true for everyone who reads it.

**A name is the shortcut and the columns are the identification.** Every block
carries a ``COLUMN_SIGNATURE``, because a file name is a convention people
depart from -- they rename what they download, and one format here has no
default name to depart from in the first place. The name is still asked first,
so a file called what its engine calls it is claimed without being opened.

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

    Some tables are too small and too plainly named to be signed by overlap.
    DIA-NN's XIC export is ``pr``, ``feature``, ``rt``, ``value`` -- four words
    that belong to nobody, and ``rt`` is JMod's as well. What tells them apart is
    not which of those names is present but that *all* of them are and there is
    almost nothing else: four columns against JMod's thirty-four.

    So both have to hold. ``columns`` are all required, which on its own already
    separates the two -- JMod has ``rt`` and none of the other three. ``width``
    is the most columns the table is expected to have, and is the second lock: a
    frame that happens to carry all four names among two hundred others is not
    this table, whatever else it is.
    """

    #: Every one of these has to be present, not two of them.
    columns: frozenset[str]
    #: The most columns the table has. Generous, since it is guarding against
    #: something an order of magnitude wider rather than measuring the table.
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
    #: too -- a frame this package renamed and wrote back out reads as DIA-NN,
    #: which is what it is. `PEP` is left out for being MaxQuant's as well.
    #:
    #: The report only. `xic` is signed by its shape instead, below.
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
    #: The XIC export, which no overlap of names could claim: `pr`, `feature`,
    #: `rt` and `value` are four words that belong to nobody and one of them is
    #: JMod's. Its shape is what identifies it -- all four present, in a table
    #: four columns wide where an identification table is thirty-four or more.
    #: The names are the ones DIA-NN is read under downstream; the ceiling is
    #: loose on purpose, being there to exclude something an order of magnitude
    #: wider rather than to measure the table.
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
    #: Fifteen tables that share few columns, so this has to reach any of them
    #: rather than describe one: `Raw file` carries most, the scan tables are
    #: signed by their instrument readings and `proteinGroups`, which has no
    #: run column at all, by its protein ones.
    #:
    #: `Charge` and `Intensity` are deliberately absent. FragPipe writes both
    #: under exactly those names, and two of them together are all it would take
    #: to make a psm.tsv ambiguous -- which is refused rather than guessed at, so
    #: it would not mis-read one file, it would make both unreadable.
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
    #: All three identification tables carry one set of columns, so one
    #: signature reaches every table this format has. `rt` and `mz` are left
    #: out: they are generic enough that DIA-NN's XIC export writes `rt` too.
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
    #: The search scores and the spectrum columns for `psm`, the counting ones
    #: for `protein`, so either table is reached. `Charge`, `Intensity`,
    #: `Peptide`, `Protein` and `Gene` are all absent for being too plain to be
    #: anybody's in particular -- the first two are MaxQuant's under the same
    #: names, and a file two blocks claim is refused rather than guessed at.
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
#: tab-separated export spells them. The parquet export writes the same report
#: with ``_`` where this has ``.``, so :class:`Spectronaut` carries both
#: spellings -- derived from this rather than written out twice, because two
#: lists of seventeen names differing by one character is two lists that drift.
#: Measured off both: an export of each was read and every key checked against
#: what the file actually holds.
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

    A rename mapping is applied by name, and a name the file does not carry
    does nothing, so holding both costs a dict twice the size and buys not
    having to know which serialization is being read.
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

    **The separator depends on the serialization.** The tab-separated export
    writes ``R.FileName``; the parquet export writes ``R_FileName``, and turns
    the space in ``PG.Cscore (Run-Wise)`` into an underscore as well. A dot is
    a path separator in a nested parquet schema, so the export spells it out of
    the way. Both are the same report and both are mapped.

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

    #: A name the export sometimes goes by, not a name it must have. Kept
    #: because a report exported as one is then recognized without opening it.
    FILES: list[str] = field(default_factory=lambda: ["Report"])
    #: ``<date>_<time>_<analysis>_Report`` is a shape an export often has and
    #: never has to: **Spectronaut has no default output name** -- whoever runs
    #: the analysis names it, and `GluC-30min.parquet` is as real an export as
    #: any. So this is a shortcut, not the identification; COLUMN_SIGNATURE is
    #: what actually settles it. Matched against the stem in full rather than
    #: from its start, so a ``..._Report.setup`` beside a report is not taken
    #: for it, and case-sensitively, because DIA-NN's ``report`` differs from a
    #: bare ``Report`` by the one letter and a file two blocks claim is an
    #: error rather than a guess.
    FILE_PATTERNS: list[str] = field(default_factory=lambda: [r".*_Report"])
    #: Both. The pattern is over the stem, so it does not care which. DIA-NN
    #: claims ``.parquet`` too, which is what the case-sensitivity above is
    #: holding apart.
    FILE_EXTENSIONS: list[str] = field(default_factory=lambda: [".parquet", ".tsv"])
    #: What identifies the report when its name cannot, which is the usual
    #: case. Any two of these settle it: they are level-prefixed the way no
    #: other engine read here prefixes anything, in either spelling, and a
    #: report configured without one of them still carries the rest. This is
    #: the same call the cellenONE reader makes for the same reason -- which
    #: file is which is worked out from the file, because names are unreliable.
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
    #: There is no EG.PrecursorId in every export -- there was none in the one
    #: this was written from -- so the identifier the rest of the package keys
    #: on is built from the two columns a fragment group always has. Stated
    #: under the canonical names, since it is built after the rename.
    BUILT_COLS: dict[str, tuple[str, ...]] = field(
        default_factory=lambda: {
            "Precursor.Id": ("Modified.Sequence", "Precursor.Charge")
        }
    )
    #: Columns the export writes as text that are not text, named as the file
    #: names them so a caller keeping the file's own names gets them too.
    #:
    #: Spectronaut is inconsistent about this within one file: ``EG.Qvalue``
    #: arrives as the string ``'1.99e-13'`` while ``PG.Qvalue`` beside it is a
    #: double, and ``PEP.IsProteotypic`` is ``'False'`` while ``EG.IsDecoy`` is
    #: a real boolean. A string q-value is not a lesser q-value, it is one that
    #: raises `TypeError` the first time anyone filters on it.
    #:
    #: A list rather than a rule, because the rule gets it wrong:
    #: ``FG.XICDBID`` is every-value-parses-as-a-number and is a database key,
    #: and turning an identifier into an integer is the kind of quiet damage
    #: this package exists not to do. ``EG.IsVerified`` is left alone for the
    #: opposite reason -- every value in the measured export is the string
    #: ``'NaN'``, so there is nothing there to say what it would be.
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
    #: As above. ``'True'`` and ``'False'``, which the tab-separated export
    #: writes as real booleans, so this is also what keeps the two agreeing.
    BOOLEAN_COLS: frozenset[str] = field(
        default_factory=lambda: _either_spelling({"PEP.IsProteotypic"})
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
