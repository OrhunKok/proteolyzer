"""The input-format descriptors must be well-formed."""

import re

from proteolyzer.core.formats import Config as CoreConfig
from proteolyzer.core.models import _claims, _shaped, _signed


def _engines(cfg) -> list[str]:
    return [name for name in vars(cfg) if hasattr(getattr(cfg, name), "FILES")]


def test_a_block_describes_the_format_and_nothing_else():
    cfg = CoreConfig()
    assert cfg.DIANN.FILES
    assert cfg.MaxQuant.COLS_RENAME_MAPPING["Experiment"] == "Run"


def test_the_format_descriptors_no_longer_nest_reference_data():
    """Settings and domain constants are separate; see proteolyzer.reference."""
    cfg = CoreConfig()
    assert not hasattr(cfg, "Protease")
    assert not hasattr(cfg, "AminoAcids")


def test_every_engine_block_is_well_formed():
    """Whatever the engines are, each describes itself the same way."""
    cfg = CoreConfig()
    engines = _engines(cfg)
    assert {"DIANN", "MaxQuant", "JMod", "FragPipe", "Spectronaut"} <= set(engines)

    for name in engines:
        block = getattr(cfg, name)
        assert block.FILES, name
        assert block.FILE_EXTENSIONS, name
        assert all(extension.startswith(".") for extension in block.FILE_EXTENSIONS), (
            name
        )
        # A pattern that does not compile would be a format nothing recognizes.
        for pattern in getattr(block, "FILE_PATTERNS", ()):
            re.compile(pattern)


def test_no_two_blocks_claim_the_same_file():
    """A pattern reaches names its author did not list, and one that reached
    another block's would make that engine's own output unreadable -- detection
    refuses to guess between two claimants. This says which pair, before a real
    file has to. DIA-NN's `report.tsv` and a bare Spectronaut `Report.tsv` are
    the near miss it is here for: one capital letter apart, same extension."""
    cfg = CoreConfig()
    engines = _engines(cfg)

    for name in engines:
        block = getattr(cfg, name)
        for stem in block.FILES:
            for extension in block.FILE_EXTENSIONS:
                claimants = [
                    other
                    for other in engines
                    if _claims(getattr(cfg, other), stem, extension)
                ]
                assert claimants == [name], (stem, extension, claimants)


def test_no_block_carries_a_column_subset():
    """Which columns to keep is the caller's, not the format's -- one list cannot
    be right for a dashboard and a pipeline at once, and carrying one here meant
    every consumer either took a subset built for somebody else or overrode it.
    ``Data.cols_to_load`` is where a project says what it wants."""
    cfg = CoreConfig()

    for name in _engines(cfg):
        assert not hasattr(getattr(cfg, name), "LOAD_COLS"), name


def test_the_maxquant_tables_name_the_run_differently():
    """evidence names the run Experiment; the rest name it Raw file."""
    cfg = CoreConfig()
    assert cfg.MaxQuant.COLS_RENAME_MAPPING["Experiment"] == "Run"
    assert cfg.MaxQuant.COLS_RENAME_MAPPING["Raw file"] == "Run"


def test_jmod_and_fragpipe_map_onto_the_canonical_names():
    cfg = CoreConfig()
    assert cfg.JMod.COLS_RENAME_MAPPING["file_name"] == "Run"
    assert cfg.JMod.COLS_RENAME_MAPPING["stripped_seq"] == "Stripped.Sequence"
    assert cfg.FragPipe.COLS_RENAME_MAPPING["Spectrum File"] == "Run"
    assert cfg.FragPipe.COLS_RENAME_MAPPING["Peptide"] == "Stripped.Sequence"


def test_spectronaut_maps_its_prefixed_names_onto_the_canonical_ones():
    """Prefixed by the level the column belongs to, which nothing else here is."""
    mapping = CoreConfig().Spectronaut.COLS_RENAME_MAPPING
    assert mapping["R.FileName"] == "Run"
    assert mapping["PEP.StrippedSequence"] == "Stripped.Sequence"
    assert mapping["FG.Charge"] == "Precursor.Charge"
    assert mapping["FG.Quantity"] == "Precursor.Quantity"
    assert mapping["EG.Qvalue"] == "Q.Value"
    assert mapping["PG.Qvalue"] == "PG.Q.Value"
    # PG.Quantity is the canonical name already, so it is not renamed at all.
    assert "PG.Quantity" not in mapping


def test_spectronaut_says_how_to_build_a_precursor_identifier():
    """Its report has no one column for it, and the rest of the package keys on
    one. Stated in the core's names, since it is built after the rename."""
    built = CoreConfig().Spectronaut.BUILT_COLS
    assert built["Precursor.Id"] == ("Modified.Sequence", "Precursor.Charge")


def test_only_spectronaut_is_recognised_by_a_pattern():
    """Every other engine names its own output; this one stamps the name with
    the date, the time and the analysis, so there is nothing fixed to list."""
    cfg = CoreConfig()
    patterned = [
        name
        for name in _engines(cfg)
        if getattr(getattr(cfg, name), "FILE_PATTERNS", ())
    ]
    assert patterned == ["Spectronaut"]


def test_spectronaut_maps_both_spellings_of_its_names():
    """The tab-separated export writes `R.FileName`, the parquet one writes
    `R_FileName`, and the block carries both -- a rename mapping is applied by
    name, so a name the file lacks does nothing and holding both costs a dict
    twice the size."""
    mapping = CoreConfig().Spectronaut.COLS_RENAME_MAPPING

    for dotted, canonical in mapping.items():
        assert mapping[dotted.replace(".", "_")] == canonical

    assert mapping["R.FileName"] == mapping["R_FileName"] == "Run"
    assert mapping["EG.Qvalue"] == mapping["EG_Qvalue"] == "Q.Value"


def test_spectronaut_can_be_identified_without_its_name():
    """It has no default output name -- the analyst names the export -- so the
    columns have to be able to say what it is on their own."""
    signature = CoreConfig().Spectronaut.COLUMN_SIGNATURE

    assert {"R.FileName", "R_FileName"} <= signature
    assert len(signature) > 2, "one column matching cannot be enough to claim a file"


def test_every_format_can_be_recognised_without_its_name():
    """A file name is a convention, not a fact about the file: people rename
    what they download. Every engine carries a signature, so a report keeps
    being read as one whatever it ends up called."""
    cfg = CoreConfig()
    for name in _engines(cfg):
        signature = getattr(getattr(cfg, name), "COLUMN_SIGNATURE", frozenset())
        assert len(signature) > 2, f"{name} cannot be recognised by its columns"


#: What each engine's tables actually hold, column for column. The DIA-NN report
#: is read off the real one in `examples/`; the rest are the per-table lists
#: `streamlit-DO-MS` keeps, which were themselves recovered from the `LOAD_COLS`
#: this package carried through v0.3.x. Written down rather than derived from the
#: signatures, because a signature checked against itself checks nothing.
REAL_TABLES: dict[str, dict[str, set[str]]] = {
    "MaxQuant": {
        "evidence": {
            "Raw file",
            "Retention time",
            "Retention length",
            "PEP",
            "Type",
            "Intensity",
            "m/z",
            "Sequence",
            "Charge",
            "MS/MS count",
            "Modified sequence",
            "PIF",
            "Missed cleavages",
            "Experiment",
            "Leading razor protein",
            "Gene names",
        },
        "allPeptides": {
            "Charge",
            "Intensity",
            "Mass",
            "Mass deficit",
            "Raw file",
            "Retention length (FWHM)",
            "Retention time",
            "Type",
            "m/z",
        },
        "msScans": {
            "Ion injection time",
            "MS/MS count",
            "Raw file",
            "Retention time",
            "Total ion current",
        },
        "msmsScans": {
            "Charge",
            "Ion injection time",
            "MS scan number",
            "Modified sequence",
            "Precursor apex offset time",
            "Raw file",
            "Sequence",
        },
        "proteinGroups": {
            "Protein IDs",
            "Majority protein IDs",
            "Potential contaminant",
            "Intensity",
        },
    },
    "JMod": {
        "filtered_IDs": {
            "BestChannel_Qvalue",
            "MS1_Area",
            "MS1_Int",
            "PredVal",
            "Protein_Qvalue",
            "Qvalue",
            "channel",
            "channels_matched",
            "file_name",
            "hyperscore",
            "is_decoy",
            "iso_cor",
            "mz",
            "pep_len",
            "plex_Area",
            "protein",
            "rt",
            "scribe_scores",
            "seq",
            "silac_channel",
            "stripped_seq",
            "tic",
            "untag_prec",
            "untag_seq",
            "window_mz",
            "z",
            "frac_dia_int",
        },
    },
    "FragPipe": {
        "psm": {
            "Assigned Modifications",
            "Calibrated Observed M/Z",
            "Charge",
            "Delta Mass",
            "Entry Name",
            "Expectation",
            "Gene",
            "Hyperscore",
            "Intensity",
            "Ion Mobility",
            "Modified Peptide",
            "Nextscore",
            "Number of Enzymatic Termini",
            "Number of Missed Cleavages",
            "Observed M/Z",
            "Peptide",
            "Peptide Length",
            "PeptideProphet Probability",
            "Protein",
            "Protein ID",
            "Retention",
            "Spectrum",
            "Spectrum File",
        },
        "protein": {
            "Coverage",
            "Description",
            "Entry Name",
            "Gene",
            "Organism",
            "Protein",
            "Protein ID",
            "Protein Probability",
            "Razor Intensity",
            "Razor Peptides",
            "Razor Spectral Count",
            "Top Peptide Probability",
            "Total Intensity",
            "Total Peptides",
            "Total Spectral Count",
            "Unique Intensity",
            "Unique Peptides",
            "Unique Spectral Count",
        },
    },
    "Spectronaut": {
        "Report": {
            "R_FileName",
            "EG_ModifiedSequence",
            "PEP_StrippedSequence",
            "FG_Charge",
            "FG_Quantity",
            "PG_ProteinGroups",
            "EG_Qvalue",
            "FG_XICDBID",
            "PG_Cscore_(Run-Wise)",
        },
    },
    "DIANN": {
        # Four words that belong to nobody, one of which is JMod's. Recognised
        # by the whole shape rather than by an overlap of names.
        "xic": {"pr", "feature", "rt", "value"},
        "report": {
            "Run.Index",
            "Run",
            "Channel",
            "Precursor.Id",
            "Modified.Sequence",
            "Stripped.Sequence",
            "Precursor.Charge",
            "Precursor.Lib.Index",
            "Proteotypic",
            "Precursor.Mz",
            "Protein.Ids",
            "Protein.Group",
            "Genes",
            "RT",
            "iRT",
            "IM",
            "Precursor.Quantity",
            "Precursor.Normalised",
            "Ms1.Area",
            "Ms1.Normalised",
            "PG.MaxLFQ",
            "Q.Value",
            "PEP",
            "Global.Q.Value",
            "Lib.Q.Value",
            "PG.Q.Value",
        },
    },
}


def _claimants(cfg, engines, columns) -> list[str]:
    """Who claims a file holding `columns`, the way `input_type` asks it: by
    whole shape first, and by an overlap of names only where nothing does."""
    return [name for name in engines if _shaped(getattr(cfg, name), columns)] or [
        name for name in engines if _signed(getattr(cfg, name), columns)
    ]


def test_every_real_table_is_claimed_by_exactly_one_format():
    """The invariant the signatures live or die by.

    Detection refuses a file two formats claim, so a signature that reached
    another engine's table would not mis-read one file -- it would make both
    unreadable. `Charge` and `Intensity` are MaxQuant's and FragPipe's under the
    same names, `PEP` is MaxQuant's and DIA-NN's, `rt` is JMod's and DIA-NN's:
    four names that cannot be in any signature, and the reason to check the
    tables rather than trust that nobody used them.
    """
    cfg = CoreConfig()
    engines = _engines(cfg)

    for engine, tables in REAL_TABLES.items():
        for table, columns in tables.items():
            claimants = _claimants(cfg, engines, columns)
            assert claimants == [engine], f"{engine}/{table} claimed by {claimants}"


def test_a_signature_is_not_satisfied_by_one_column_alone():
    """One familiar name turns up in frames people derive from a report and
    write back out; two together are what an engine's own output has."""
    cfg = CoreConfig()
    for name in _engines(cfg):
        for column in getattr(cfg, name).COLUMN_SIGNATURE:
            assert not _signed(getattr(cfg, name), {column}), column


def test_a_narrow_table_needs_every_one_of_its_columns():
    """Two of four would be an overlap, and an overlap is what cannot tell this
    table from anybody else's -- `rt` alone is JMod's too."""
    cfg = CoreConfig()
    xic = REAL_TABLES["DIANN"]["xic"]

    assert _shaped(cfg.DIANN, xic)
    for missing in xic:
        assert not _shaped(cfg.DIANN, xic - {missing}), missing


def test_a_narrow_table_is_narrow_or_it_is_not_that_table():
    """The second lock. A frame carrying all four of those words among two
    hundred others is something else, whatever else it is."""
    cfg = CoreConfig()
    xic = REAL_TABLES["DIANN"]["xic"]

    assert not _shaped(cfg.DIANN, xic | {f"extra{n}" for n in range(20)})
