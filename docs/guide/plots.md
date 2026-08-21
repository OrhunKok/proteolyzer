# Plotting

Publication-styled plots built on seaborn and matplotlib. Importing
`proteolyzer` does not import matplotlib; the module loads on first access.

```python
import proteolyzer as pz

plot = pz.plots.VolcanoPlot(
    results,
    x="log2FC",
    y="pvalue",
    label="Gene",
    signif=0.05,
    effect_threshold=1.0,
)
plot.label_points()
plot.add_data_point_count()
plot.save("volcano.png")
```

`VolcanoPlot` classifies each point as `up`, `down` or `notsig` from the
significance threshold and, if given, a minimum absolute effect size — a tiny
fold change with a great p-value is not a hit. Keyword arguments other than the
documented ones are forwarded to `seaborn.scatterplot`.

## The y axis

`y` may hold raw p-values or `-log10` p-values. Which one is inferred from the
data, but the heuristic cannot tell an untransformed column whose best p-value
is 0.5 from a transformed one whose best is 10^-0.5, so say so explicitly when
it matters:

```python
pz.plots.VolcanoPlot(results, x="log2FC", y="neglog10p", transformed=True)
```

A p-value of exactly 0 — from a permutation test, or underflow — would be
`inf` after the transform, which matplotlib drops from the axes while still
counting the point as significant: tallied in the box but invisible in the
plot. Those are drawn at the most extreme finite value instead, and the
substitution is logged.

## Themes

The default `science` theme renders text with LaTeX and so needs a LaTeX
installation. Without one, pass another theme:

```python
pz.plots.VolcanoPlot(results, x="log2FC", y="pvalue", theme="default")
```

An unknown theme is reported and ignored rather than raising, so a missing
stylesheet does not lose you the figure.
