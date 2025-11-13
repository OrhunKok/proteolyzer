---
sidebar_label: matrix
title: proteolyzer.transformers.matrix
---

## pd

## np

## Callable

## MetaLogging

## ProcessedData

## MatrixBuilder Objects

```python
class MatrixBuilder(metaclass=MetaLogging)
```

#### \_\_slots\_\_

#### \_\_init\_\_

```python
def __init__(processed_data: ProcessedData)
```

#### \_missingness\_check

```python
def _missingness_check(matrix: pd.DataFrame,
                       warning_threshold: float = 0.75) -> None
```

#### matrix\_generation

```python
def matrix_generation(values: str, index: list[str],
                      columns: list[str]) -> pd.DataFrame
```

#### normalize\_matrix

```python
def normalize_matrix(within_groups: list[str],
                     agg_func: Callable,
                     replace_zeros: bool = True) -> pd.DataFrame
```

