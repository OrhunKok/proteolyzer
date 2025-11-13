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

## np

## sns

## pd

## Literal

## List

## Rectangle

## adjust\_text

## PlotBase

## RelPlot Objects

```python
class RelPlot(PlotBase)
```

#### \_\_init\_\_

```python
def __init__(**kwargs)
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
             hue_order: List = None,
             label: str = None,
             signif: float = 0.05,
             **kwargs)
```

#### \_prepare\_data

```python
def _prepare_data() -> pd.DataFrame
```

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

