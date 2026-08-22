"""Reading a cellenONE run directory: cells, droplets, labels and the chamber.

The instrument writes a directory per run and the interesting part of it is
implicit: which log belongs to which step of the preparation, which table is the
cells and which the droplets it rejected, and which of two runs of the same step
is the one that counted. :class:`CoordinatesMapping` works that out from the
files themselves rather than from their names.
"""

from .mapping import (
    CELL_KEYS,
    CELLEONE_MAPPING,
    DISPENSING_STAGES,
    STAGES,
    UNREAD,
    CoordinatesMapping,
)

__all__ = [
    "CELLEONE_MAPPING",
    "CELL_KEYS",
    "DISPENSING_STAGES",
    "STAGES",
    "UNREAD",
    "CoordinatesMapping",
]
