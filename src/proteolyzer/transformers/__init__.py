"""Data transformation helpers for proteolyzer.

This package groups matrix and tensor transformation utilities used for
normalization, filtering, and reshaping of proteomics data.

Submodules
    matrix: matrix-specific transforms (normalization, scaling, imputation)
"""

from .matrix import MatrixBuilder

__all__ = ["MatrixBuilder"]
