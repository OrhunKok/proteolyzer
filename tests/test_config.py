"""The input-format descriptors must be well-formed."""

import re

from proteolyzer.core.formats import Config as CoreConfig
from proteolyzer.core.models import _claims


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
