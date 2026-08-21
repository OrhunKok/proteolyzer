"""Proteolyzer: processing, analysis and visualization of proteomics data.

The core (``config``, ``utils``, ``transformers``) is imported eagerly. The
optional subpackages listed in ``__lazy__`` are imported on first attribute
access so that a missing extra (``pip install proteolyzer[aas]``) only fails
for the module that needs it.
"""

import importlib

from . import config, transformers, utils
from .utils import Data, configure_logging

__all__ = [
    "Data",
    "aas",
    "cellenone",
    "config",
    "configure_logging",
    "plots",
    "transformers",
    "unimod",
    "utils",
]

#: Subpackages imported on first attribute access, see :func:`__getattr__`.
__lazy__ = ("aas", "cellenone", "plots", "unimod")


def __getattr__(name: str):
    if name in __lazy__:
        module = importlib.import_module(f"{__name__}.{name}")
        globals()[name] = module
        return module

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__lazy__))
