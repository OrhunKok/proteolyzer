"""The shared reference data must be complete, lazy and immutable."""

import pytest

from proteolyzer import reference


def test_amino_acid_masses():
    masses = reference.amino_acid_masses()
    assert masses["A"] == pytest.approx(71.03711, abs=1e-4)
    assert masses["G"] == pytest.approx(57.02146, abs=1e-4)


def test_three_letter_codes():
    codes = reference.three_letter_to_one()
    assert codes["Ala"] == "A"
    # Xle is added by the UniMod export step for isobaric Ile/Leu.
    assert codes["Xle"] == "J"


def test_the_modification_table_is_available():
    mods = reference.modifications()
    assert {"one_letter", "mono_mass", "classification", "full_name"} <= set(
        mods.columns
    )
    assert "AA substitution" in set(mods["classification"])


def test_tables_are_read_once():
    assert reference.amino_acids() is reference.amino_acids()
    assert reference.amino_acid_masses() is reference.amino_acid_masses()


def test_standard_codon_table_is_complete():
    assert len(reference.CODON_TABLE) == 64
    assert reference.CODON_TABLE["AUG"] == "M"
    assert set(reference.CODON_TABLE.values()) >= set("ACDEFGHIKLMNPQRSTVWY*")


@pytest.mark.parametrize(
    ("name", "sites", "counts"),
    [
        ("Trypsin", ("K", "R"), {"K": 1, "R": 1}),
        ("LysC", ("K",), {"K": 1}),
        ("ArgC", ("R",), {"R": 1}),
    ],
)
def test_proteases(name, sites, counts):
    protease = reference.protease(name)
    assert protease.name == name
    assert protease.cleavage_sites == sites
    assert dict(protease.allowed_counts) == counts


def test_an_unknown_protease_lists_the_valid_ones():
    with pytest.raises(
        ValueError, match=r"Must be one of: \['ArgC', 'LysC', 'Trypsin'\]"
    ):
        reference.protease("Chymotrypsin")


def test_reference_data_cannot_be_mutated():
    """Regression: `sites += [...]` extended the shared list on every call."""
    with pytest.raises(TypeError):
        reference.protease("Trypsin").allowed_counts["K"] = 9
    with pytest.raises(TypeError):
        reference.PROTEASES["Trypsin"] = None
    with pytest.raises(TypeError):
        reference.CODON_TABLE["AUG"] = "X"
    with pytest.raises(TypeError):
        reference.amino_acid_masses()["A"] = 0.0

    assert reference.protease("Trypsin").cleavage_sites == ("K", "R")


def test_a_missing_resource_is_reported():
    with pytest.raises(FileNotFoundError, match="nope.csv"):
        reference.resource_stream("nope.csv")
