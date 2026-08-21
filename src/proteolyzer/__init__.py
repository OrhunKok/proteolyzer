"""Proteolyzer: processing, analysis and visualization of proteomics data.

The base suite is :mod:`proteolyzer.core` — recognizing search-engine output,
reading it, and normalizing it into a consistent frame — with
:mod:`proteolyzer.reference` holding the domain constants they share. Both are
imported eagerly.

Everything else is an optional subpackage, imported on first attribute access
so that a missing extra (``pip install proteolyzer[unimod]``) only fails for
the module that needs it. ``tests/test_package_boundaries.py`` keeps that true.
"""

import importlib

from . import core, reference
from .core import (
    Data,
    DataLoader,
    DataProcessor,
    MatrixBuilder,
    Report,
    configure_logging,
    read,
)

__all__ = [
    "Data",
    "DataLoader",
    "DataProcessor",
    "MatrixBuilder",
    "Report",
    "configure_logging",
    "core",
    "plots",
    "read",
    "reference",
    "unimod",
]

#: Subpackages imported on first attribute access, see :func:`__getattr__`.
__lazy__ = ("plots", "unimod")


def __getattr__(name: str):
    if name in __lazy__:
        module = importlib.import_module(f"{__name__}.{name}")
        globals()[name] = module
        return module

    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__lazy__))
