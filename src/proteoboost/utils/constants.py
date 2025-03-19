# DIA-NN related constants
DIANN_FILES = [
    'report',
    'report.stats',
    'report.log',
]
DIANN_FILES += [f"{file}-first-pass" for file in DIANN_FILES]

DIANN_COLS = {
    'report.stats': None,  # Load all columns
    'report': {
        'Run', 'Precursor.Id', 'Stripped.Sequence', 'Precursor.Charge', 'Genes',
        'Protein.Group', 'Proteotypic', 'RT', 'RT.Start', 'RT.Stop', 'FWHM',
        'IM', 'PEP', 'Ms1.Area', 'Ms1.Apex.Area', 'Ms1.Normalised',
        'Precursor.Quantity', 'Precursor.Normalised',
    },
    'xic': None,  # Load all columns
}

DIANN_EXTENSIONS = ['.parquet', '.tsv']

# MaxQuant related constants
MAXQUANT_FILES = [
    'allPeptides', 'evidence', 'matchedFeatures', 'modificationSpecificPeptides.txt',
    'ms3Scans', 'msms', 'msmsScans', 'mzRange', 'Oxidation (M)Sites',
    'parameters', 'peptides', 'proteinGroups', 'summary',
]

MAXQUANT_COLS = {
    'allPeptides': None,  # Load all columns
    'evidence': {
        'Experiment', 'Modified sequence', 'Sequence', 'Charge', 'Gene names',
        'Proteins', 'Leading razor protein', 'Missed cleavages',
        'Retention time', 'Retention length', 'PEP', 'Intensity',
    },
    'msScans': None,  # Load all columns
    'msmsScans': None,  # Load all columns
    'parameters': None,  # Load all columns
    'summary': None,  # Load all columns
}

MAXQUANT_EXTENSIONS = ['.txt']

# Supported files and columns
SUPPORTED_FILES = {'DIANN': DIANN_FILES, 'MaxQuant': MAXQUANT_FILES}
SUPPORTED_FILES_COLS_SUBSET = {'DIANN': DIANN_COLS, 'MaxQuant': MAXQUANT_COLS}

# Column renaming mappings
COLS_RENAME_MAPPING = {
    'DIANN': {},
    'MaxQuant': {
        'Experiment': 'Run',
        'Modified sequence': 'Precursor.Id',
        'Sequence': 'Stripped.Sequence',
        'Charge': 'Precursor.Charge',
        'Gene names': 'Genes',
        'Protein.Group': 'Protein.Group',
        'Leading razor protein': 'Leading.Razor.Protein',
        'Missed cleavages': 'Missed.Cleavages',
        'Retention time': 'RT',
        'Retention length': 'RT.Width',
    },
}

# Processing constants
COL_MEDIAN_THRESHOLD = 100
CARDINALITY_THRESHOLD = 0.1

# Protease rules
PROTEASE_RULES = {
    'Trypsin': {'K': 1, 'R': 1},
    'LysC': {'K': 1},
    'ArgC': {'R': 1},
}