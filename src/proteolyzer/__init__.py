"""Top-level package for proteolyzer.

This package collects subpackages that handle data loading, transformations,
visualizations and utilities used throughout the proteolyzer project.

Docstring style
    Google-style docstrings are used across the project to maximize
    compatibility with pydoc-markdown and other documentation tools.

Package contents
    - cellenone: I/O and processing specific to CellenONE data
    - plots: visualization helpers and plotting abstractions
    - transformers: matrix and tensor transformations
    - utils: constants, helpers, and core processors
"""

from . import utils
from . import plots
from . import transformers

__all__ = ["utils", "plots", "transformers"]
