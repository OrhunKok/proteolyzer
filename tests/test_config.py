"""The input-format descriptors must be well-formed."""

from proteolyzer.core.formats import Config as CoreConfig


def test_core_format_blocks_are_dicts_and_sets():
    cfg = CoreConfig()
    assert isinstance(cfg.DIANN.LOAD_COLS, dict)
    assert isinstance(cfg.DIANN.LOAD_COLS["report"], set)
    assert cfg.DIANN.LOAD_COLS["report.stats"] is None
    assert cfg.DIANN.LOAD_COLS["report"] == cfg.DIANN.LOAD_COLS["report-first-pass"]
    assert cfg.MaxQuant.COLS_RENAME_MAPPING["Experiment"] == "Run"


def test_the_format_descriptors_no_longer_nest_reference_data():
    """Settings and domain constants are separate; see proteolyzer.reference."""
    cfg = CoreConfig()
    assert not hasattr(cfg, "Protease")
    assert not hasattr(cfg, "AminoAcids")


def test_every_engine_block_is_well_formed():
    """Whatever the engines are, each describes itself the same way."""
    cfg = CoreConfig()
    engines = [name for name in vars(cfg) if hasattr(getattr(cfg, name), "FILES")]
    assert {"DIANN", "MaxQuant", "JMod", "FragPipe"} <= set(engines)

    for name in engines:
        block = getattr(cfg, name)
        assert block.FILES, name
        assert block.FILE_EXTENSIONS, name
        assert set(block.LOAD_COLS) <= set(block.FILES), name
        assert all(
            columns is None or isinstance(columns, set)
            for columns in block.LOAD_COLS.values()
        ), name
        assert all(extension.startswith(".") for extension in block.FILE_EXTENSIONS), (
            name
        )


def test_the_maxquant_tables_besides_evidence_have_subsets():
    """A file with no subset is read whole, and allPeptides runs to several GB."""
    cfg = CoreConfig()
    for name in ("allPeptides", "msScans", "msmsScans"):
        assert isinstance(cfg.MaxQuant.LOAD_COLS[name], set), name
        assert "Raw file" in cfg.MaxQuant.LOAD_COLS[name], name

    # evidence names the run Experiment; the rest name it Raw file.
    assert cfg.MaxQuant.COLS_RENAME_MAPPING["Raw file"] == "Run"


def test_jmod_and_fragpipe_map_onto_the_canonical_names():
    cfg = CoreConfig()
    assert cfg.JMod.COLS_RENAME_MAPPING["file_name"] == "Run"
    assert cfg.JMod.COLS_RENAME_MAPPING["stripped_seq"] == "Stripped.Sequence"
    assert cfg.FragPipe.COLS_RENAME_MAPPING["Spectrum File"] == "Run"
    assert cfg.FragPipe.COLS_RENAME_MAPPING["Peptide"] == "Stripped.Sequence"

    # Both write their identifications under one name per file, so a subset is
    # an intersection with whatever that workflow wrote.
    assert "plex_Area" in cfg.JMod.LOAD_COLS["filtered_IDs"]
    assert "Hyperscore" in cfg.FragPipe.LOAD_COLS["psm"]
