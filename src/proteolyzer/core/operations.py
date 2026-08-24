"""Small pure operations used throughout proteolyzer.

This module holds focused, well-documented pure functions that operate on
core in-memory data representations (lists, dicts, DataFrames).
"""

from collections.abc import Callable, Sequence

import numpy as np
import pandas as pd


def per_distinct(
    func: Callable[[pd.Series], pd.Series],
) -> Callable[[pd.Series], pd.Series]:
    """Wrap a column transformation so it only runs on the distinct values.

    ``func`` must be *element-wise*: the result for a row may depend on that
    row's value and nothing else. Given that, applying it to the distinct
    values and gathering the results back out is the same answer for far less
    work -- pandas runs string operations element by element in Python, even
    on categorical columns, and a report repeats every peptide, protein group
    and identifier once per run per channel.

    Missing values are passed through ``func`` rather than short-circuited, so
    however it treats them is preserved. The one difference from applying
    ``func`` row by row: factorizing normalizes the missing sentinel, so a
    ``None`` in a column of nothing but ``None`` comes back as ``NaN``. The
    dtype is unchanged and both read as missing.

    >>> s = pd.Series(["B;A", "B;A", "C"])
    >>> per_distinct(lambda x: x.str.split(";").str[0])(s).tolist()
    ['B', 'B', 'C']
    """

    def wrapper(column: pd.Series) -> pd.Series:
        # use_na_sentinel=False keeps NA as a value of its own rather than as
        # code -1, so func decides what it means and the result dtype is not
        # widened to hold a fill value that never gets used.
        codes, uniques = pd.factorize(column, use_na_sentinel=False)
        derived = func(pd.Series(pd.Index(uniques)))
        gathered = derived.reindex(codes)
        gathered.index = column.index
        return gathered

    return wrapper


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


def log10(values: Sequence[float] | np.ndarray) -> np.ndarray:
    """log10 of `values`, with zeros and negative numbers dropped first.

    Written for an axis: ``np.log10`` of zero is ``-inf`` and of a negative
    number is ``nan``, either of which sends a plot's scale to a value with
    nothing to show for it. Dropping the offending values first, rather than
    clipping them, keeps what survives exact.

    >>> log10([1.0, 10.0, 0.0, -5.0]).tolist()
    [0.0, 1.0]
    """
    data = np.asarray(values, dtype="float64")
    return np.log10(data[data > 0])


def jaccard_index(
    left: Sequence[bool] | np.ndarray, right: Sequence[bool] | np.ndarray
) -> float:
    """How much two masks agree: what they share over what either of them has.

    Written for asking how far two labelling channels, or two runs, identified
    the same precursors -- where the answer wanted is the overlap rather than
    either count on its own.

    Both are read as boolean masks over the same ordered set, so they have to be
    the same length and aligned. Returns ``nan`` for two empty masks, where the
    ratio is 0/0 and there is nothing to be similar about.

    >>> jaccard_index([True, True, False], [True, False, False])
    0.5
    """
    left = np.asarray(left, dtype=bool)
    right = np.asarray(right, dtype=bool)
    union = int(np.logical_or(left, right).sum())
    if union == 0:
        return np.nan
    return float(np.logical_and(left, right).sum() / union)
