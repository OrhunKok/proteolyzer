from dataclasses import dataclass
from typing import ClassVar


CODON_MAP = {
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

TMT_PLEX_MAP = {'126' : '1', '127N' : '2', '127C' : '3', '128N' : '4', '128C' : '5', '129N' : '6',
                '129C' : '7', '130N' : '8', '130C' : '9', '131N' : '10', '131C' : '11', '132N' : '12',
                '132C' : '13', '133N' : '14', '133C' : '15', '134N' : '16', '134C' : '17', '135N' : '18'}

MQ_TMT_MAP = {'126' : 1,'127N' : 2, '127C' : 3, '128N' : 4, 
               '128C' : 5, '129N' : 6, '129C' : 7, '130N' : 8, 
               '130C' : 9, '131' : 10}

PROTEASE_CLEAVAGE_SITES = {'Trypsin' : ['K', 'R', '*'], 'Lys-C' : ['K', '*'], 'Arg-C' : ['R', '*']}


@dataclass(frozen=True)
class AAS:
    
    TMT_PLEX_MAP: ClassVar[dict[str, str]] = {
        '126' : '1', '127N' : '2', '127C' : '3', '128N' : '4', '128C' : '5', 
        '129N' : '6', '129C' : '7', '130N' : '8', '130C' : '9', '131N' : '10', 
        '131C' : '11', '132N' : '12', '132C' : '13', '133N' : '14', '133C' : '15', 
        '134N' : '16', '134C' : '17', '135N' : '18'
    }

    MQ_TMT_MAP: ClassVar[dict[str, int]] = {
        '126' : 1, '127N' : 2, '127C' : 3, '128N' : 4, 
        '128C' : 5, '129N' : 6, '129C' : 7, '130N' : 8, 
        '130C' : 9, '131' : 10
    }
