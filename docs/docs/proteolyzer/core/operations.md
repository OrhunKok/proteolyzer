---
sidebar_label: operations
title: proteolyzer.core.operations
---

Small pure operations used throughout proteolyzer.

This module holds focused, well-documented pure functions that operate on
core in-memory data representations (lists, dicts, DataFrames).

## Callable

## Sequence

## np

## pd

#### per\_distinct

```python
def per_distinct(
    func: Callable[[pd.Series],
                   pd.Series]) -> Callable[[pd.Series], pd.Series]
```

Wrap a column transformation so it only runs on the distinct values.

``func`` must be *element-wise*: the result for a row may depend on that
row&#x27;s value and nothing else. Given that, applying it to the distinct
values and gathering the results back out is the same answer for far less
work -- pandas runs string operations element by element in Python, even
on categorical columns, and a report repeats every peptide, protein group
and identifier once per run per channel.

Missing values are passed through ``func`` rather than short-circuited, so
however it treats them is preserved. The one difference from applying
``func`` row by row: factorizing normalizes the missing sentinel, so a
``None`` in a column of nothing but ``None`` comes back as ``NaN``. The
dtype is unchanged and both read as missing.

&gt;&gt;&gt; s = pd.Series([&quot;B;A&quot;, &quot;B;A&quot;, &quot;C&quot;])
&gt;&gt;&gt; per_distinct(lambda x: x.str.split(&quot;;&quot;).str[0])(s).tolist()
[&#x27;B&#x27;, &#x27;B&#x27;, &#x27;C&#x27;]

#### cv

```python
def cv(data: Sequence[float] | np.ndarray, min_datapoints: int = 3) -> float
```

Coefficient of variation (sample standard deviation / mean).

Returns ``nan`` when there are fewer than ``min_datapoints`` observations or
when the mean is zero, so the result can be aggregated without guarding.

