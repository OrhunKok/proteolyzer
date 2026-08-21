"""Amino acid substitution (AAS) discovery pipeline.

Implements the pipeline used for discovery of amino acid substitutions and
PTMs (https://decode.slavovlab.net/). The stages are meant to be run in order,
each reading the same parameter file:

    Preprocessor -> FrameTranslator -> Detection -> Validation -> Quantification

Modules
    base: shared stage plumbing (parameter loading, queue, metadata)
    params: parameter file schema and loader
    results: read back what a run produced
    preprocessing: search-engine output to parquet conversion
    translation: six-frame translation of a genome FASTA
    detection: SAAP and ALT assignment
    validation: fragment-level validation of candidates
    quantification: substitution ratio quantification
"""

from .base import NullQueue, Stage
from .detection import Detection, MaxQuant
from .params import ParamsSchema, load_params
from .preprocessing import Preprocessor
from .quantification import Quantification
from .results import ARTEFACTS, Results
from .translation import FrameTranslator
from .validation import Validation

__all__ = [
    "ARTEFACTS",
    "Detection",
    "FrameTranslator",
    "MaxQuant",
    "NullQueue",
    "ParamsSchema",
    "Preprocessor",
    "Quantification",
    "Results",
    "Stage",
    "Validation",
    "load_params",
]
