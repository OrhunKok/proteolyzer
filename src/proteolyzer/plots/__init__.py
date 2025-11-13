"""Plotting utilities and factories for proteolyzer.

This package exposes plotting helpers and base classes used to render
diagnostic and relational plots for proteomics datasets.

Submodules
    base: common plotting base classes and utilities
    relational: relational plot implementations (pairwise, scatter, etc.)

"""

from .relational import VolcanoPlot

__all__ = ["VolcanoPlot"]
