---
sidebar_label: io
title: proteolyzer.core.io
---

Interchange format for frames passed between pipeline stages.

Stages used to hand each other pickled DataFrames. Pickle ties a results
folder to the pandas version that produced it and cannot be read safely from
elsewhere, which is a poor property for output archived alongside a paper, so
frames are written as parquet instead.

This lives in the core rather than in one subpackage: persisting a frame is
not specific to any pipeline, and a subpackage should not have to import
another one to do it.

## Path

## pd

#### FRAME\_SUFFIX

#### frame\_path

```python
def frame_path(path: Path) -> Path
```

The parquet path for `path`.

``Path.with_suffix`` replaces everything after the last dot, so a sample
named ``pt.01_SAAP`` would become ``pt.parquet`` — and two samples could
collapse onto the same file. Only this module&#x27;s own suffix is replaced.

#### write\_frame

```python
def write_frame(df: pd.DataFrame, path: Path) -> Path
```

Write `df` as parquet, returning the path actually written.

#### read\_frame

```python
def read_frame(path: Path) -> pd.DataFrame
```

Read a frame written by :func:`write_frame`.

#### frame\_exists

```python
def frame_exists(path: Path) -> bool
```

Whether a frame is available at `path`.

