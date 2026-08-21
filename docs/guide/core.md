# Loading and processing

## Reading a report

```python
import proteolyzer as pz

report = pz.read("report.parquet")
```

`read` recognizes known DIA-NN and MaxQuant file names and reads only the
columns configured for them. Use `load_all_columns=True` for everything, or
`extra_cols_to_load=[...]` to add to the default set.

```python
pz.read("report.tsv", load_all_columns=True)
pz.read("report.tsv", extra_cols_to_load=["Q.Value", "Predicted.RT"])
```

Parquet, TSV, CSV, Excel and plaintext are dispatched on the file extension,
and a file-like object works in place of a path as long as it has a `name` to
dispatch on. Which search engine produced a file is detected from its name and
extension, and can be overridden with `INPUT_TYPE=`.

Delimited files are read with pyarrow's multithreaded parser where its peak
memory is affordable, and with the stock parser otherwise — see
[Memory and dtypes](memory.md#reading).

## The Report

[`Report`][proteolyzer.core.models.Report] is a frozen wrapper around three
things: the `frame`, the `source` it was read from, and the `processing` that
produced it (`None` until `process()` runs).

```python
report.frame                    # the pandas DataFrame
report.source.input_type        # 'DIANN'
report.columns
len(report)
report["Run"]
```

It is deliberately **not** a `DataFrame` subclass. pandas returns plain frames
from most operations, so a subclass loses its metadata the first time anyone
slices it, and the internals it would have to hook into are not a stable API.
`len()`, `report["Run"]` and `.columns` are there for interactive use;
everything else goes through `.frame`, which is always a plain `DataFrame`.

Being frozen, a report's metadata cannot drift from the frame it describes. To
carry the metadata onto a different frame, ask for it:

```python
subset = report.with_frame(report.frame.query("Q.Value < 0.01"))
```

## Processing

```python
processed = report.process()
```

`process()` normalizes the frame and records what it did:

- drops columns holding the same value in every row
- narrows dtypes to what the values need (see [Memory and dtypes](memory.md))
- renames columns to the canonical proteolyzer names
- derives `Leading.Razor.Protein`, `Peptide.Length`, `RT.Width` and a
  label-free identifier
- detects labelled precursors and derives channel information from the
  identifiers
- flags peptides whose residue counts are inconsistent with full cleavage

The input report is untouched; a new one comes back.

```python
processed.processing.label_free        # what processing found
processed.processing.labels_complete   # False if channels could not be derived
processed.n_identifications            # distinct precursors
```

### Labelling

Labelled data is detected by matching the precursor identifier against a
regex, which defaults to mTRAQ, SILAC and TMT. For each label found, processing
adds `<label>.Label`, `<label>.Count` and `<label>.Channel`, plus
`Run.<label>.Channel` and `Run.Full.Channel`.

If a channel cannot be assigned — offsets that are not uniform across a
peptide, say — the columns are left out and `processing.labels_complete` is
`False`. Check it before relying on channel-level results.

```python
processed.process(label_group_capture=r"\((mTRAQ[^()]*)\)")   # narrow it
```

### Identification counts

```python
processed.summary()
```

|  | Rows | Precursors | Peptides | Proteins |
|---|---|---|---|---|
| dOK083 | 315 | 315 | 315 | 286 |
| dOK084 | 300 | 300 | 298 | 272 |

Counts are of *distinct* values, so a precursor seen in two channels of one run
counts once. Levels the frame does not carry are left out.

## Matrices

```python
import numpy as np

builder = processed.matrix("Ms1.Area", index=["Precursor.Id"], columns=["Run"])
matrix = builder.normalize_matrix(within_groups=["Run"], agg_func=np.nansum).matrix
```

`normalize_matrix` divides each value by an aggregate of its row within each
column group, which is how a plex is normalized within a run.

Missingness is measured when the matrix is built, and the figures are returned
rather than only logged:

```python
missing = builder.missingness()

missing.mar                  # % of cells that are NA
missing.mnar                 # % that are exactly zero
missing.per_column           # fraction absent per column
missing.sparse_columns(0.75) # candidates for dropping
```

Gaps and zeros are counted apart because they mean different things: an NA is a
measurement that was never made, a zero is one that came back empty, and
imputation has to treat them differently.
