"""Small pure operations used throughout proteolyzer.

This module holds focused, well-documented pure functions that operate on
core in-memory data representations (lists, dicts, DataFrames).
"""

from collections.abc import Sequence

import numpy as np


def cv(data: Sequence[float] | np.ndarray, min_datapoints: int = 3) -> float:
    """Coefficient of variation (sample standard deviation / mean).

    Returns ``nan`` when there are fewer than ``min_datapoints`` observations or
    when the mean is zero, so the result can be aggregated without guarding.
    """
    data = np.asarray(data, dtype="float64")
    if data.size < min_datapoints:
        return np.nan
    mean = np.mean(data)
    if mean == 0:
        return np.nan
    return np.std(data, ddof=1) / mean
