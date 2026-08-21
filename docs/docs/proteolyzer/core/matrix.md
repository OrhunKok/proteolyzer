---
sidebar_label: matrix
title: proteolyzer.core.matrix
---

Matrix transformation utilities.

Provide performant, well-documented routines to perform common matrix
operations required by the analysis pipeline, including normalization,
scaling and basic imputations.

Example
    &gt;&gt;&gt; import numpy as np
    &gt;&gt;&gt; builder = MatrixBuilder(processed_data)  # doctest: +SKIP
    &gt;&gt;&gt; builder.matrix_generation(
    ...     &quot;Ms1.Area&quot;, index=[&quot;Precursor.Id&quot;], columns=[&quot;Run&quot;]
    ... ).normalize_matrix(within_groups=[&quot;Run&quot;], agg_func=np.nansum).matrix

## Callable

## dataclass

## np

## pd

## Logged

## Report

## Missingness Objects

```python
@dataclass(frozen=True)
class Missingness()
```

How much of a matrix is absent, and from where.

Gaps and zeros are counted apart because they mean different things: an NA
is a measurement that was never made (missing at random), a zero is one
that came back empty (missing not at random), and imputation has to treat
them differently.

#### mar

Percentage of cells that are NA.

#### mnar

Percentage of cells that are exactly zero.

#### per\_column

Fraction of each column that is absent, counting both gaps and zeros.

#### sparse\_columns

```python
def sparse_columns(threshold: float = 0.75) -> list
```

Columns more than `threshold` absent, as candidates for dropping.

## MatrixBuilder Objects

```python
class MatrixBuilder(Logged)
```

Pivots long-form processed data into a quantitative matrix.

#### \_\_slots\_\_

#### matrix

Set by :meth:`matrix_generation`.

#### \_\_init\_\_

```python
def __init__(data: Report | pd.DataFrame)
```

Takes a :class:`Report` or a plain frame.

#### missingness

```python
def missingness(matrix: pd.DataFrame | None = None,
                warning_threshold: float = 0.75) -> Missingness
```

Measure how much of the matrix is absent, logging the headline numbers.

Defaults to the generated matrix. Returns the numbers as well as
logging them, so a caller can act on them rather than read them.

#### matrix\_generation

```python
def matrix_generation(values: str, index: list[str],
                      columns: list[str]) -> MatrixBuilder
```

#### \_require\_matrix

```python
def _require_matrix() -> pd.DataFrame
```

#### normalize\_matrix

```python
def normalize_matrix(within_groups: list[str],
                     agg_func: Callable,
                     replace_zeros: bool = True) -> MatrixBuilder
```

Divide each value by an aggregate of its row within each column group.

Parameters
----------
within_groups
    Column-index level(s) defining the groups aggregated over.
agg_func
    Row-wise aggregation, called as ``agg_func(block, axis=1)`` (e.g.
    ``np.nansum``, ``np.nanmedian``).
replace_zeros
    Treat zeros as missing before normalizing.

