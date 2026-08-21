---
sidebar_label: operations
title: proteolyzer.core.operations
---

Small pure operations used throughout proteolyzer.

This module holds focused, well-documented pure functions that operate on
core in-memory data representations (lists, dicts, DataFrames).

## Sequence

## np

#### cv

```python
def cv(data: Sequence[float] | np.ndarray, min_datapoints: int = 3) -> float
```

Coefficient of variation (sample standard deviation / mean).

Returns ``nan`` when there are fewer than ``min_datapoints`` observations or
when the mean is zero, so the result can be aggregated without guarding.

