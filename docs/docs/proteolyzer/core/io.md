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

:func:`read_frame` still accepts the old ``.p`` files, so a half-finished run
from a previous version can be picked up where it left off.

## pickle

## Path

## pd

#### FRAME\_SUFFIX

#### LEGACY\_SUFFIX

#### \_with\_suffix

```python
def _with_suffix(path: Path, suffix: str) -> Path
```

`path` carrying `suffix`, without eating part of the name.

``Path.with_suffix`` replaces everything after the last dot, so a sample
named ``pt.01_MTP`` would become ``pt.parquet`` — and two samples could
collapse onto the same file. Only a suffix this module owns is replaced.

#### frame\_path

```python
def frame_path(path: Path) -> Path
```

The parquet path for `path`, whatever suffix it was given.

#### legacy\_path

```python
def legacy_path(path: Path) -> Path
```

The pre-parquet pickle path for `path`.

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

