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
