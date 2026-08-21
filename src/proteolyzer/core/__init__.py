"""The base suite: read search-engine output into a consistent frame.

Submodules
    formats: the input formats recognized, and how their columns are handled
    loader: read a described source into memory
    models: Data, LoadedData and ProcessedData
    processor: dtype narrowing, derived columns, labelling information
    matrix: pivot processed data into a quantitative matrix
    operations: small pure functions
    logging: the package logger and the Logged base class
    io: parquet interchange for frames passed between stages
    pipeline: shared stage plumbing (parameters, progress, provenance)
"""

from .io import frame_exists, read_frame, write_frame
from .loader import DataLoader
from .logging import Logged, configure_logging
from .matrix import MatrixBuilder
from .models import Data, LoadedData, ProcessedData
from .operations import cv
from .pipeline import NullQueue, Stage
from .processor import DataProcessor

__all__ = [
    "Data",
    "DataLoader",
    "DataProcessor",
    "LoadedData",
    "Logged",
    "MatrixBuilder",
    "NullQueue",
    "ProcessedData",
    "Stage",
    "configure_logging",
    "cv",
    "frame_exists",
    "read_frame",
    "write_frame",
]
