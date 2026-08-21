# Proteolyzer

![Proteolyzer](assets/proteolyzer_logo.svg){ width="420" }

A Python package for processing, analyzing and visualizing proteomics data. It
reads DIA-NN and MaxQuant output, normalizes it into a consistent shape, and
provides the domain pipelines used for single-cell sample preparation and amino
acid substitution discovery.

## What it does

- **Data loading** — Parquet, TSV/CSV, Excel and plaintext, with per-format
  column subsetting so only the columns you need are read, through pyarrow's
  multithreaded parser where it can be used.
- **Data processing** — dtype narrowing, derived columns, missed-cleavage
  flags, per-run identification counts, and automatic detection of labelled
  (mTRAQ/SILAC/TMT) precursors.
- **Matrix transformation** — pivot to a quantitative matrix, measure
  missingness, and normalize within column groups.
- **Plotting** — publication-styled relational plots.
- **[Reference data](guide/reference-data.md)** — amino acid masses, the
  genetic code and protease rules, plus SQL access to the full UniMod database,
  which is built on first use rather than shipped.

## Quickstart

```python
import numpy as np
import proteolyzer as pz

report = pz.read("report.parquet")   # DIA-NN report, recognized by name
processed = report.process()          # dtypes, derived columns, labelling info

processed.summary()                   # precursors/peptides/proteins per run
processed.memory()                    # what each column costs
processed.frame                       # the pandas DataFrame

matrix = (
    processed.matrix("Ms1.Area", index=["Precursor.Id"], columns=["Run"])
    .normalize_matrix(within_groups=["Run"], agg_func=np.nansum)
    .matrix
)
```

Start with the [core guide](guide/core.md), or go straight to the
[API reference](reference/proteolyzer/index.md).

## Design

The core — `proteolyzer.core` and `proteolyzer.reference` — is imported with
the package and depends only on pandas, numpy and pyarrow. `plots` and
`unimod` load on first access, so a missing optional dependency only fails for
the module that needs it. `tests/test_package_boundaries.py` enforces that,
including that importing proteolyzer does not pull in matplotlib.

Domain pipelines built on this core live in their own repositories, so an
instrument or an assay does not become a dependency of everyone's install:

- [proteolyzer-cellenone](https://github.com/OrhunKok/proteolyzer-cellenone) —
  mapping single cells prepared on a cellenONE to well positions
- [proteolyzer-aas](https://github.com/OrhunKok/proteolyzer-aas) — discovery of
  amino acid substitutions
