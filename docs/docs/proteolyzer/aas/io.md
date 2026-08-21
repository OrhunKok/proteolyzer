---
sidebar_label: io
title: proteolyzer.aas.io
---

Interchange format for the frames passed between AAS pipeline stages.

Stages used to hand each other pickled DataFrames. Pickle ties a results
folder to the pandas version that produced it and cannot be read safely from
elsewhere, which is a poor property for output archived alongside a paper, so
frames are written as parquet instead.

:func:`read_frame` still accepts the old ``.p`` files, so a half-finished run
from a previous version can be picked up where it left off.

## pickle

## Path

## pd

#### FRAME\_SUFFIX

#### LEGACY\_SUFFIX

#### frame\_path

```python
def frame_path(path: Path) -> Path
```

The parquet path for `path`, whatever suffix it was given.

#### write\_frame

```python
def write_frame(df: pd.DataFrame, path: Path) -> Path
```

Write `df` as parquet, returning the path actually written.

#### read\_frame

```python
def read_frame(path: Path) -> pd.DataFrame
```

Read a frame written by :func:`write_frame`, or a legacy pickle.

#### frame\_exists

```python
def frame_exists(path: Path) -> bool
```

Whether a frame is available at `path`, in either format.

