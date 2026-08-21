"""Amino acid substitution (AAS) discovery pipeline.

Implements the pipeline used for discovery of amino acid substitutions
(SAAP) and the modifications that would otherwise explain them (ALT); see
https://decode.slavovlab.net/. The stages are meant to be run in order,
each reading the same parameter file:

    Preprocessor -> FrameTranslator -> Detection -> Validation -> Quantification

Modules
    base: shared stage plumbing (parameter loading, queue, metadata)
    pipeline: run the stages in order, in the two phases the method allows
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
from .pipeline import Pipeline
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
    "Pipeline",
    "Preprocessor",
    "Quantification",
    "Results",
    "Stage",
    "Validation",
    "load_params",
]
