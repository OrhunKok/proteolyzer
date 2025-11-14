MAXQUANT_FILES = []
DIANN_FILES = []

PEP_THRESHOLD = 0.01

AA_MW = {
        "G": 57.02147, "A" : 71.03712, "S" : 87.03203, "P" : 97.05277, "V" : 99.06842, "T" : 101.04768,
        "I" : 113.08407, "L" : 113.08407, "N" : 114.04293, "D" : 115.02695, "Q" : 128.05858, "K" : 128.09497,
        "E" : 129.0426,"M" : 131.04049, "H" : 137.05891, "F" : 147.06842, "R" : 156.10112, "C" : 160.030654, 
        "Y" : 163.0633, "W" : 186.07932,
        }
        
TMT_MASS_SHIFT = {0 : 224.152478, 2 : 225.155833, 6 : 229.162932, 10 : 229.162932, 11: 229.169252}


COLSTOKEEP = {
            'MaxQuant': {
                        'allpeptides.txt' : [
                                'Raw file', 'Charge', 'm/z', 'Mass', 'Retention time',
                                'Intensity', 'DP Mass Difference', 'DP PEP', 'DP Decoy', 'DP Proteins',
                                'Reverse', 'DP Base Sequence', 'DP Probabilities', 'DP Positional Probability',
                                'DP Base Scan Number', 'DP Mod Scan Number', 'MSMS Scan Numbers'
                            ],
                        'evidence.txt' : [
                                            'Raw file', 'Charge', 'm/z', 'Mass', 'Retention time', 'Reverse', 'Potential contaminant',
                                            'Sequence', 'PIF', 'PEP', 'Mass error [ppm]', 'MS/MS scan number', 'Intensity'
                                        ],
                        'msms.txt' : ['Raw file', 'Scan number', 'Matches', 'Reverse'],
                        'peptides.txt' : ['Sequence', 'Start position', 'Amino acid after', 'Amino acid before',
                                 'Reverse', 'Potential contaminant']
                         }
            }

FILES_NEEDED = {
            'MaxQuant': {'Detection' : ['evidence.txt', 'allPeptides.txt', 'peptides.txt'], 'Validation' : ['evidence.txt', 'msms.txt']}
            }