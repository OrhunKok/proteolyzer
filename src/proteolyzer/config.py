from dataclasses import dataclass
from typing import Dict, Final, ClassVar, List
import pandas as pd
import importlib.resources
import io
from pathlib import Path


_RESOURCE_MODULE: Final[str] = 'proteolyzer.resources'
_AA_FILE: Final[str] = "unimod_amino_acids.csv"

def _get_resource_stream(filename: str) -> io.BytesIO:
    resource_path: Path = importlib.resources.files(_RESOURCE_MODULE).joinpath(filename)
    try:
        file_bytes = resource_path.read_bytes()
    except FileNotFoundError as e:
        raise FileNotFoundError(f"Resource file '{filename}' not found.") from e
    
    return io.BytesIO(file_bytes)


@dataclass(frozen=True)
class AminoAcids:
    """
    Container for amino acid properties.
    """
    MASS: ClassVar[Final[Dict[str, float]]] = (
        pd.read_csv(
            _get_resource_stream(_AA_FILE),
            usecols=["one_letter", "mono_mass"],
        )
        .set_index("one_letter")["mono_mass"]
        .to_dict()
    )
    CODE_TO_SYMBOL: ClassVar[Final[Dict[str, str]]] = (
        pd.read_csv(
            _get_resource_stream(_AA_FILE),
            usecols=["three_letter", "one_letter"],
        )
        .set_index("three_letter")["one_letter"]
        .to_dict()
    )


@dataclass(frozen=True)
class Codons:
    """
    Container for genetic code mappings, organized by species.
    """
    Standard: ClassVar[Dict[str, str]] = {
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
        "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"
    }



@dataclass(frozen=True)
class Protease:
    """
    Defines cleavage rules and allowed residue counts for common proteases.
    Each protease is a nested subclass with:
        - CLEAVAGE_SITES: list of residues where cleavage occurs
        - ALLOWED_COUNTS: dictionary representing maximum allowed counts of specific residues in peptides.
    """

    @dataclass(frozen=True)
    class Trypsin:
        CLEAVAGE_SITES: ClassVar[List[str]] = ["K", "R"]
        ALLOWED_COUNTS: ClassVar[Dict[str, int]] = {"K": 1, "R": 1}
        

    @dataclass(frozen=True)
    class LysC:
        CLEAVAGE_SITES: ClassVar[List[str]] = ["K"]
        ALLOWED_COUNTS: ClassVar[Dict[str, int]] = {"K": 1}

    @dataclass(frozen=True)
    class ArgC:
        CLEAVAGE_SITES: ClassVar[List[str]] = ["R"]
        ALLOWED_COUNTS: ClassVar[Dict[str, int]] = {"R": 1}