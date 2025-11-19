"""Top-level package for proteolyzer.

This package collects subpackages that handle data loading, transformations,
visualizations and utilities used throughout the proteolyzer project.

Package contents
    - plots: visualization helpers and plotting abstractions
    - transformers: matrix and tensor transformations
    - utils: constants, helpers, and core processors
"""

from . import config
from . import utils
from . import plots
from . import transformers
from . import plugins

from .utils import Data


__all__ = ["config", "utils", "plots", "transformers", "plugins", 'Data']
