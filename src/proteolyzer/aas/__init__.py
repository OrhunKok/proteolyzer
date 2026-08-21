"""Amino acid substitution (AAS) discovery pipeline.

Implements the pipeline used for discovery of amino acid substitutions and
PTMs (https://decode.slavovlab.net/). The stages are meant to be run in order,
each reading the same parameter file:

    Preprocessor -> FrameTranslator -> Detection -> Validation -> Quantification

Modules
    base: shared stage plumbing (parameter loading, queue, metadata)
    params: parameter file schema and loader
    preprocessing: search-engine output to parquet conversion
    translation: six-frame translation of a genome FASTA
    detection: candidate substitution and PTM detection
    validation: fragment-level validation of candidates
    quantification: substitution ratio quantification
"""

from .base import NullQueue, Stage
from .detection import Detection, MaxQuant
from .params import ParamsSchema, load_params
from .preprocessing import Preprocessor
from .quantification import Quantification
from .translation import FrameTranslator
from .validation import Validation

__all__ = [
    "Detection",
    "FrameTranslator",
    "MaxQuant",
    "NullQueue",
    "ParamsSchema",
    "Preprocessor",
    "Quantification",
    "Stage",
    "Validation",
    "load_params",
]
