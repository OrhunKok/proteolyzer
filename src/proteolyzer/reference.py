"""Domain reference data: one source of truth for the constants.

Everything here is fixed knowledge about proteins and proteases rather than a
setting, so it is exposed as module-level immutable mappings and cached
loaders instead of a configuration object. Subpackages import what they need
directly; they do not re-expose it through their own settings objects, which
only added a layer to look through.

The tables come from the UniMod export in ``proteolyzer/resources``, read on
first use rather than at import, so importing the package costs nothing. See
``proteolyzer.unimod`` for how those files are regenerated.
"""

import functools
import importlib.resources
import io
from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType
from typing import Final

import pandas as pd

_RESOURCE_MODULE: Final[str] = "proteolyzer.resources"
AMINO_ACID_FILE: Final[str] = "unimod_amino_acids.csv"
MODIFICATION_FILE: Final[str] = "unimod_modifications.csv"


def resource_stream(filename: str) -> io.BytesIO:
    """The bundled resource `filename`, as a stream."""
    resource = importlib.resources.files(_RESOURCE_MODULE).joinpath(filename)
    try:
        return io.BytesIO(resource.read_bytes())
    except FileNotFoundError as e:
        raise FileNotFoundError(f"Resource file '{filename}' not found.") from e


@functools.cache
def amino_acids() -> pd.DataFrame:
    """The amino acid table, including Xle (``J``) for isobaric Ile/Leu."""
    return pd.read_csv(resource_stream(AMINO_ACID_FILE))


@functools.cache
def modifications() -> pd.DataFrame:
    """Every approved UniMod modification, one row per residue it applies to."""
    return pd.read_csv(resource_stream(MODIFICATION_FILE))


@functools.cache
def amino_acid_masses() -> Mapping[str, float]:
    """Monoisotopic residue mass by one-letter code."""
    table = amino_acids().set_index("one_letter")["mono_mass"]
    return _frozen(table.to_dict())


@functools.cache
def three_letter_to_one() -> Mapping[str, str]:
    """One-letter code by three-letter code, e.g. ``{"Ala": "A"}``."""
    table = amino_acids().set_index("three_letter")["one_letter"]
    return _frozen(table.to_dict())


def _frozen(mapping: dict) -> Mapping:
    """A read-only view, so a shared table cannot be mutated by a caller."""
    return MappingProxyType(mapping)


#: The standard genetic code, keyed by RNA codon. ``*`` is a stop codon.
CODON_TABLE: Final[Mapping[str, str]] = _frozen({
    "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
    "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
})  # fmt: skip


@dataclass(frozen=True)
class Protease:
    """Cleavage rules for one protease.

    Attributes
    ----------
    cleavage_sites
        Residues the protease cleaves after.
    allowed_counts
        How many times each of those residues may appear in a fully cleaved
        peptide. A tuple and a read-only mapping, because these are shared:
        extending them in place used to corrupt every later call.
    """

    name: str
    cleavage_sites: tuple[str, ...]
    allowed_counts: Mapping[str, int]


#: Proteases the processing pipeline can flag missed cleavages for.
PROTEASES: Final[Mapping[str, Protease]] = _frozen(
    {
        protease.name: protease
        for protease in (
            Protease("Trypsin", ("K", "R"), _frozen({"K": 1, "R": 1})),
            Protease("LysC", ("K",), _frozen({"K": 1})),
            Protease("ArgC", ("R",), _frozen({"R": 1})),
        )
    }
)


def protease(name: str) -> Protease:
    """The named protease, or a ValueError listing the ones there are."""
    try:
        return PROTEASES[name]
    except KeyError:
        raise ValueError(
            f"Invalid protease: '{name}'. Must be one of: {sorted(PROTEASES)}"
        ) from None
