from dataclasses import dataclass, field
from typing import ClassVar, List, Dict
from proteolyzer import config as CONFIG


@dataclass(frozen=True) 
class DIANN:
    FILES: List[str] = field(default_factory=lambda: [
        "report", "report-first-pass"
    ])
    LOAD_COLS: Dict[str, set] = field(default_factory=lambda: {
        "report.stats": None,  # Load all columns
        "report": None, # Load all columns
        "report-first-pass": None # Load all columnss
    })
    FILE_EXTENSIONS: List[str] = field(default_factory=lambda: [".parquet", ".tsv"])
    COLS_RENAME_MAPPING: Dict[str, str] = field(default_factory=lambda: {})
    FILES_NEEDED: Dict[str, set] = field(default_factory=lambda: {
        "Detection": ['report'],
        "Validation": ['report'],
    })


@dataclass(frozen=True)
class MaxQuant:
    FILES: List[str] = field(default_factory=lambda: [
        "allPeptides", "evidence", "msms", "peptides"
    ])
    LOAD_COLS: Dict[str, set] = field(default_factory=lambda: {
        "allPeptides": [
                        'Raw file', 'Charge', 'm/z', 'Mass', 'Retention time',
                        'Intensity', 'DP Mass Difference', 'DP PEP', 'DP Decoy', 'DP Proteins',
                        'Reverse', 'DP Base Sequence', 'DP Probabilities', 'DP Positional Probability',
                        'DP Base Scan Number', 'DP Mod Scan Number', 'MSMS Scan Numbers'
                        ],
        "evidence": [
                'Raw file', 'Charge', 'm/z', 'Mass', 'Retention time', 'Reverse', 'Potential contaminant',
                'Sequence', 'PIF', 'PEP', 'Mass error [ppm]', 'MS/MS scan number', 'Intensity'
                ],
        "msms": ['Raw file', 'Scan number', 'Matches', 'Reverse'],
        "peptides": ['Sequence', 'Start position', 'Amino acid after', 'Amino acid before',
                                 'Reverse', 'Potential contaminant'],
    })
    FILE_EXTENSIONS: List[str] = field(default_factory=lambda: [".txt"])
    COLS_RENAME_MAPPING: Dict[str, str] = field(default_factory=lambda: {
        "Experiment": "Run",
        "Modified sequence": "Precursor.Id",
        "Sequence": "Stripped.Sequence",
        "Charge": "Precursor.Charge",
        "Gene names": "Genes",
        "Protein.Group": "Protein.Group",
        "Leading razor protein": "Leading.Razor.Protein",
        "Missed cleavages": "Missed.Cleavages",
        "Retention time": "RT",
        "Retention length": "RT.Width",
    }),
    FILES_NEEDED: Dict[str, set] = field(default_factory=lambda: {
        "Detection": ['evidence', 'allPeptides', 'peptides'],
        "Validation": ['evidence', 'msms'],
    })


@dataclass(frozen=True)
class TMT:
    
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
    MASS_SHIFT: ClassVar[dict[int, float]] = {0 : 224.152478, 2 : 225.155833, 6 : 229.162932, 10 : 229.162932, 11: 229.169252}



@dataclass(frozen=True)
class Config:
    PEP_THRESHOLD: float = 0.01
    DIANN: DIANN = field(default_factory=DIANN)
    MaxQuant: MaxQuant = field(default_factory=MaxQuant)
    TMT: TMT = field(default_factory=TMT)
    Protease: CONFIG.Protease = field(default_factory=CONFIG.Protease)
    AminoAcids: CONFIG.Protease = field(default_factory=CONFIG.AminoAcids)
