"""Utility helpers and common constants for proteolyzer.

This package contains low-level utilities that are reused across the
project, such as configuration/constants, data loaders, logging helpers,
and small processing utilities.

Submodules
    - config: recognized input formats and their column handling
    - loader: functions to load common file types
    - logging: logging configuration and helpers
    - models: small data classes and typed models
    - operations: small pure functions operating on core data structures
    - processor: higher-level processing pipelines
"""

from .loader import DataLoader
from .logging import configure_logging
from .models import Data, LoadedData, ProcessedData
from .operations import cv
from .processor import DataProcessor

__all__ = [
    "Data",
    "DataLoader",
    "DataProcessor",
    "LoadedData",
    "ProcessedData",
    "configure_logging",
    "cv",
]
