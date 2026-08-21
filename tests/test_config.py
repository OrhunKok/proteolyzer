"""The input-format descriptors must be well-formed."""

from proteolyzer.aas.config import Config as AASConfig
from proteolyzer.core.formats import Config as CoreConfig


def test_core_format_blocks_are_dicts_and_sets():
    cfg = CoreConfig()
    assert isinstance(cfg.DIANN.LOAD_COLS, dict)
    assert isinstance(cfg.DIANN.LOAD_COLS["report"], set)
    assert cfg.DIANN.LOAD_COLS["report.stats"] is None
    assert cfg.DIANN.LOAD_COLS["report"] == cfg.DIANN.LOAD_COLS["report-first-pass"]
    assert cfg.MaxQuant.COLS_RENAME_MAPPING["Experiment"] == "Run"


def test_aas_config_field_types():
    """Regression: a stray trailing comma made COLS_RENAME_MAPPING a tuple."""
    cfg = AASConfig()
    assert isinstance(cfg.MaxQuant.COLS_RENAME_MAPPING, dict)
    assert cfg.MaxQuant.COLS_RENAME_MAPPING["Sequence"] == "Stripped.Sequence"
    assert isinstance(cfg.MaxQuant.FILES_NEEDED, dict)
    assert cfg.MaxQuant.FILES_NEEDED["Detection"] == [
        "evidence",
        "allPeptides",
        "peptides",
    ]


def test_the_format_descriptors_no_longer_nest_reference_data():
    """Settings and domain constants are separate; see proteolyzer.reference."""
    for cfg in (CoreConfig(), AASConfig()):
        assert not hasattr(cfg, "Protease")
        assert not hasattr(cfg, "AminoAcids")


def test_tmt_channel_map():
    assert AASConfig().TMT.MQ_TMT_MAP["126"] == 1
