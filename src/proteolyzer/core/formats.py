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
    #: Columns to load per file. ``None`` loads every column.
    LOAD_COLS: dict[str, set | None] = field(
        default_factory=lambda: {
            #: What a QC dashboard reads out of the three tables below. They had
            #: no subset, and a file with no subset is read whole -- an
            #: allPeptides.txt runs to several GB. From streamlit-DO-MS, which
            #: reads all four of these; see #23.
            "allPeptides": {
                "Raw file",
                "Retention length (FWHM)",
                "Intensity",
                "Charge",
                "m/z",
                "Retention time",
                "Type",
                "Mass",
                "Mass deficit",
            },
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
            "msScans": {
                "Raw file",
                "Ion injection time",
                "Total ion current",
                "Retention time",
                "MS/MS count",
            },
            "msmsScans": {
                "Raw file",
                "Precursor apex offset time",
                "Ion injection time",
                "Sequence",
                "Charge",
                "MS scan number",
                "Modified sequence",
            },
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


#: Columns loaded from a JMod identification table. The parquet carries the
#: first block and the csv every column, so a subset is an intersection: a plot
#: that needs one of the csv-only columns checks it arrived.
_JMOD_ID_COLS: frozenset[str] = frozenset(
    {
        "file_name",
        "seq",
        "stripped_seq",
        "z",
        "untag_prec",
        "untag_seq",
        "channel",
        "silac_channel",
        "protein",
        "is_decoy",
        "Qvalue",
        "Protein_Qvalue",
        "BestChannel_Qvalue",
        "PredVal",
        "plex_Area",
        "rt",
        "mz",
        "coeff",
        "MS1_Int",
        "MS1_Area",
        "ms1_cor",
        "iso_cor",
        "tic",
        "window_mz",
        "mz_error",
        "rt_error",
        "med_frag_error",
        "num_lib",
        "frac_lib_int",
        "frac_dia_int",
        "hyperscore",
        "scribe_scores",
        "cosine",
        "spec_r2",
        "prec_r2",
        "n_scans",
        "smoothness",
        "channels_matched",
        "frac_shared_intensity",
        "pep_len",
        "b_counts",
        "y_counts",
        "mean_frag_corr",
        "median_frag_corr",
        "frac_corr_above_0p5",
        "fitted_spectral_contrasts",
        "gof_stats",
    }
)


@dataclass(frozen=True)
class JMod:
    FILES: list[str] = field(
        default_factory=lambda: ["filtered_IDs", "all_IDs", "all_IDs_filtered"]
    )
    #: all_IDs is everything the search considered, decoys included, which is
    #: what a score distribution is read from; filtered_IDs is what survived.
    LOAD_COLS: dict[str, set | None] = field(
        default_factory=lambda: {
            "filtered_IDs": set(_JMOD_ID_COLS),
            "all_IDs": set(_JMOD_ID_COLS),
            "all_IDs_filtered": set(_JMOD_ID_COLS),
        }
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


#: Columns loaded from a FragPipe/Philosopher psm.tsv. A FragPipe table has a
#: column for every tool that ran and no two workflows run the same set, so this
#: is an intersection too: no IonQuant means no Intensity, no timsTOF means no
#: Ion Mobility. 'Probability' is what Philosopher called PeptideProphet
#: Probability before FragPipe renamed it.
_FRAGPIPE_PSM_COLS: frozenset[str] = frozenset(
    {
        "Spectrum",
        "Spectrum File",
        "Peptide",
        "Modified Peptide",
        "Peptide Length",
        "Charge",
        "Retention",
        "Observed Mass",
        "Calibrated Observed Mass",
        "Observed M/Z",
        "Calibrated Observed M/Z",
        "Calculated Peptide Mass",
        "Calculated M/Z",
        "Delta Mass",
        "Expectation",
        "Hyperscore",
        "Nextscore",
        "PeptideProphet Probability",
        "Probability",
        "Number of Enzymatic Termini",
        "Number of Missed Cleavages",
        "Intensity",
        "Ion Mobility",
        "Assigned Modifications",
        "Purity",
        "Is Unique",
        "Protein",
        "Protein ID",
        "Entry Name",
        "Gene",
        "Protein Description",
    }
)

_FRAGPIPE_PROTEIN_COLS: frozenset[str] = frozenset(
    {
        "Protein",
        "Protein ID",
        "Entry Name",
        "Gene",
        "Description",
        "Protein Probability",
        "Top Peptide Probability",
        "Total Peptides",
        "Unique Peptides",
        "Razor Peptides",
        "Total Spectral Count",
        "Unique Spectral Count",
        "Razor Spectral Count",
        "Total Intensity",
        "Unique Intensity",
        "Razor Intensity",
        "Coverage",
        "Protein Length",
        "Organism",
    }
)


@dataclass(frozen=True)
class FragPipe:
    FILES: list[str] = field(
        default_factory=lambda: ["psm", "peptide", "ion", "protein"]
    )
    LOAD_COLS: dict[str, set | None] = field(
        default_factory=lambda: {
            "psm": set(_FRAGPIPE_PSM_COLS),
            "protein": set(_FRAGPIPE_PROTEIN_COLS),
            "peptide": None,
            "ion": None,
        }
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
