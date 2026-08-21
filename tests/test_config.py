"""Configuration objects must be well-formed and free of shared mutable state."""

import pytest

from proteolyzer import config
from proteolyzer.aas.config import Config as AASConfig
from proteolyzer.utils.config import Config as CoreConfig


def test_amino_acid_tables_loaded_from_resources():
    assert config.AminoAcids.MASS["A"] == pytest.approx(71.03711, abs=1e-4)
    assert config.AminoAcids.CODE_TO_SYMBOL["Ala"] == "A"
    # Xle is added by the UniMod export step.
    assert config.AminoAcids.CODE_TO_SYMBOL["Xle"] == "J"


def test_standard_codon_table_is_complete():
    assert len(config.Codons.Standard) == 64
    assert config.Codons.Standard["AUG"] == "M"
    assert set(config.Codons.Standard.values()) >= set("ACDEFGHIKLMNPQRSTVWY*")


@pytest.mark.parametrize(
    ("protease", "sites"),
    [("Trypsin", ("K", "R")), ("LysC", ("K",)), ("ArgC", ("R",))],
)
def test_protease_cleavage_sites(protease, sites):
    assert getattr(config.Protease, protease).CLEAVAGE_SITES == sites


def test_cleavage_sites_cannot_be_extended_in_place():
    """Regression: `sites += [...]` used to mutate shared class state."""
    sites = config.Protease.Trypsin.CLEAVAGE_SITES
    with pytest.raises(TypeError):
        sites += ["*"]
    assert config.Protease.Trypsin.CLEAVAGE_SITES == ("K", "R")


def test_core_config_blocks_are_dicts_and_sets():
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


def test_aas_config_exposes_shared_reference_data():
    cfg = AASConfig()
    assert cfg.AminoAcids.MASS["G"] == pytest.approx(57.02146, abs=1e-4)
    assert cfg.Protease.Trypsin.ALLOWED_COUNTS == {"K": 1, "R": 1}
    assert cfg.TMT.MQ_TMT_MAP["126"] == 1
