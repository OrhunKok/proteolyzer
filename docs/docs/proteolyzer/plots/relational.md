---
sidebar_label: relational
title: proteolyzer.plots.relational
---

Relational plotting utilities (scatter, pairwise, joint plots).

This module implements higher-level relational plots used for exploratory
data analysis in proteolyzer. Plots accept data frames or arrays and
return Matplotlib axes or Seaborn objects.

Example
    &gt;&gt;&gt; from proteolyzer.plots.relational import scatter_plot
    &gt;&gt;&gt; ax = scatter_plot(df, x=&quot;intensity&quot;, y=&quot;ratio&quot;)

## Literal

## np

## pd

## sns

## adjust\_text

## Rectangle

## PlotBase

## RelPlot Objects

```python
class RelPlot(PlotBase)
```

Shared annotation helpers for relational (x/y) plots.

#### \_\_init\_\_

```python
def __init__(theme: str = "science", **kwargs)
```

#### \_symmetric\_xaxis

```python
def _symmetric_xaxis() -> None
```

#### label\_points

```python
def label_points(signif_filter: bool = True,
                 size: float = 4,
                 max_label: int = 50,
                 ha="center",
                 va="center",
                 **kwargs) -> None
```

#### add\_data\_point\_count

```python
def add_data_point_count(x_pos: float = 0.7625,
                         y_pos: float = 0.825,
                         fontsize: int = 6,
                         ha: str = "left",
                         va: str = "top",
                         color="black",
                         **kwargs)
```

## VolcanoPlot Objects

```python
class VolcanoPlot(RelPlot)
```

#### \_\_init\_\_

```python
def __init__(data: pd.DataFrame,
             x: str,
             y: str,
             hue: Literal["Regulation", "Significance", str] = "Regulation",
             hue_order: list = None,
             label: str = None,
             signif: float = 0.05,
             effect_threshold: float = 0.0,
             transformed: bool | None = None,
             symmetric_x: bool = False,
             delta_text_size: int = 6,
             theme: str = "science",
             **kwargs)
```

Parameters
----------
signif
    Significance threshold, as a p-value.
effect_threshold
    Minimum absolute value of `x` for a point to count as regulated.
    Zero, the default, classifies on significance alone.
transformed
    Whether `y` already holds -log10 p-values. Inferred from the data
    when None, which cannot tell an untransformed column whose best
    p-value is 0.5 from a transformed one whose best is 10^-0.5.
symmetric_x
    Centre the x axis on zero.

Remaining keyword arguments go to ``seaborn.scatterplot``. The options
above are named rather than read out of them, because anything left in
kwargs is forwarded to the plotting call.

#### \_prepare\_data

```python
def _prepare_data() -> pd.DataFrame
```

#### \_minus\_log10

```python
def _minus_log10(p_values: pd.Series) -> pd.Series
```

-log10 of `p_values`, with exact zeros held just off the top.

A p-value of exactly 0 (permutation tests, or underflow) gives inf,
which matplotlib drops from the axes while still counting the point as
significant: it would be tallied in the box but invisible in the plot.

#### \_add\_threshold\_lines

```python
def _add_threshold_lines() -> None
```

#### \_add\_delta\_count\_box

```python
def _add_delta_count_box(box_position=(0.675, 0.85),
                         box_width=0.3,
                         box_height=0.125) -> None
```

