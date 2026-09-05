"""The base suite: read search-engine output into a consistent frame.

Submodules
    formats: the input formats recognized, and how their columns are handled
    loader: read a described source into memory
    models: Data, Report and the processing metadata
    processor: dtype narrowing, derived columns, labelling information
    matrix: pivot processed data into a quantitative matrix
    isolation: what a DIA window design did to a precursor's envelope
    operations: small pure functions
    logging: the package logger and the Logged base class
    io: parquet interchange for frames passed between stages
    pipeline: shared stage plumbing (parameters, progress, provenance)
"""

from .io import frame_exists, read_frame, write_frame
from .isolation import envelope_room, envelope_split
from .loader import DataLoader
from .logging import Logged, configure_logging
from .matrix import MatrixBuilder
from .models import Data, Processing, Report
from .operations import cv, jaccard_index
from .pipeline import NullQueue, Stage
from .processor import DataProcessor, Narrower, narrow
from .reader import read

__all__ = [
    "Data",
    "DataLoader",
    "DataProcessor",
    "Logged",
    "MatrixBuilder",
    "Narrower",
    "NullQueue",
    "Processing",
    "Report",
    "Stage",
    "configure_logging",
    "cv",
    "envelope_room",
    "envelope_split",
    "frame_exists",
    "jaccard_index",
    "narrow",
    "read",
    "read_frame",
    "write_frame",
]
