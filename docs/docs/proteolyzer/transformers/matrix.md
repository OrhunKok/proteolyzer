---
sidebar_label: matrix
title: proteolyzer.transformers.matrix
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

## np

## pd

## Logged

## ProcessedData

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
def __init__(processed_data: ProcessedData)
```

#### missingness\_check

```python
def missingness_check(matrix: pd.DataFrame,
                      warning_threshold: float = 0.75) -> None
```

#### matrix\_generation

```python
def matrix_generation(values: str, index: list[str],
                      columns: list[str]) -> "MatrixBuilder"
```

#### normalize\_matrix

```python
def normalize_matrix(within_groups: list[str],
                     agg_func: Callable,
                     replace_zeros: bool = True) -> "MatrixBuilder"
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
