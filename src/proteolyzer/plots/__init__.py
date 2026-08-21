"""Plotting utilities and factories for proteolyzer.

This package exposes plotting helpers and base classes used to render
diagnostic and relational plots for proteomics datasets.

Submodules
    base: common plotting base classes and utilities
    relational: relational plot implementations (pairwise, scatter, etc.)

"""

from .base import PlotBase
from .relational import RelPlot, VolcanoPlot

__all__ = ["PlotBase", "RelPlot", "VolcanoPlot"]
