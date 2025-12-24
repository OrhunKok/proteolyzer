"""CellenONE-related helpers and entry points.

This subpackage contains utilities and wrappers to load, validate and
process CellenONE export formats used by proteolyzer.

Modules
    cellenone: high-level data loader and adapter for CellenONE exports.

"""

from .cellenone import CELLEONE_MAPPING, CoordinatesMapping

__all__ = ["CELLEONE_MAPPING", "CoordinatesMapping"]
