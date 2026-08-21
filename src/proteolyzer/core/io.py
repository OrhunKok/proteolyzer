"""Interchange format for frames passed between pipeline stages.

Stages used to hand each other pickled DataFrames. Pickle ties a results
folder to the pandas version that produced it and cannot be read safely from
elsewhere, which is a poor property for output archived alongside a paper, so
frames are written as parquet instead.

This lives in the core rather than in one subpackage: persisting a frame is
not specific to any pipeline, and a subpackage should not have to import
another one to do it.
"""

from pathlib import Path

import pandas as pd

FRAME_SUFFIX = ".parquet"


def frame_path(path: Path) -> Path:
    """The parquet path for `path`.

    ``Path.with_suffix`` replaces everything after the last dot, so a sample
    named ``pt.01_SAAP`` would become ``pt.parquet`` — and two samples could
    collapse onto the same file. Only this module's own suffix is replaced.
    """
    path = Path(path)
    if path.suffix == FRAME_SUFFIX:
        return path
    return path.with_name(path.name + FRAME_SUFFIX)


def write_frame(df: pd.DataFrame, path: Path) -> Path:
    """Write `df` as parquet, returning the path actually written."""
    target = frame_path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    # pyarrow is a core dependency and handles the list-valued columns the
    # detection and validation stages produce.
    df.to_parquet(target, engine="pyarrow")
    return target


def read_frame(path: Path) -> pd.DataFrame:
    """Read a frame written by :func:`write_frame`."""
    target = frame_path(path)
    if not target.exists():
        raise FileNotFoundError(f"No frame at {target}")
    # The default engine resolves to pyarrow, which is a core dependency;
    # naming it explicitly trips a pandas-stubs overload.
    return pd.read_parquet(target)


def frame_exists(path: Path) -> bool:
    """Whether a frame is available at `path`."""
    return frame_path(path).exists()
