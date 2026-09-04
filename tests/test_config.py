"""The input-format descriptors must be well-formed."""

import re

from proteolyzer.core.formats import Config as CoreConfig


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
    """Whatever the engines are, each describes itself the same way.

    FILES is how most of them say which files they claim; a block with no
    exact name to give -- Spectronaut's carries a timestamp -- says so with
    FILE_PATTERN instead, so a block must have one or the other.
    """
    cfg = CoreConfig()
    engines = [name for name in vars(cfg) if hasattr(getattr(cfg, name), "FILES")]
    assert {"DIANN", "MaxQuant", "JMod", "FragPipe", "Spectronaut"} <= set(engines)

    for name in engines:
        block = getattr(cfg, name)
        assert block.FILES or getattr(block, "FILE_PATTERN", ""), name
        assert block.FILE_EXTENSIONS, name
        assert all(extension.startswith(".") for extension in block.FILE_EXTENSIONS), (
            name
        )


def test_no_block_carries_a_column_subset():
    """Which columns to keep is the caller's, not the format's -- one list cannot
    be right for a dashboard and a pipeline at once, and carrying one here meant
    every consumer either took a subset built for somebody else or overrode it.
    ``Data.cols_to_load`` is where a project says what it wants."""
    cfg = CoreConfig()
    engines = [name for name in vars(cfg) if hasattr(getattr(cfg, name), "FILES")]

    for name in engines:
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


def test_spectronaut_maps_onto_the_canonical_names():
    cfg = CoreConfig()
    assert cfg.Spectronaut.COLS_RENAME_MAPPING["R.FileName"] == "Run"
    assert cfg.Spectronaut.COLS_RENAME_MAPPING["PEP.StrippedSequence"] == (
        "Stripped.Sequence"
    )
    assert cfg.Spectronaut.COLS_RENAME_MAPPING["FG.Charge"] == "Precursor.Charge"
    # Precursor.Id has no source column here -- EG.ModifiedSequence and
    # FG.Charge would have to be combined to build one -- so it is not a
    # rename and is left out.
    assert "Precursor.Id" not in cfg.Spectronaut.COLS_RENAME_MAPPING.values()


def test_spectronaut_file_pattern_matches_a_timestamped_report():
    """No exact name can match every export, so it is a pattern, not a name."""
    cfg = CoreConfig()
    assert not cfg.Spectronaut.FILES
    assert re.fullmatch(
        cfg.Spectronaut.FILE_PATTERN,
        "20260901_164751_2026-08-27_CF_PD_GluC_30min_Report",
    )
    assert not re.fullmatch(cfg.Spectronaut.FILE_PATTERN, "report")
