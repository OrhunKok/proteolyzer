"""Interchange format for frames passed between pipeline stages.

Stages used to hand each other pickled DataFrames. Pickle ties a results
folder to the pandas version that produced it and cannot be read safely from
elsewhere, which is a poor property for output archived alongside a paper, so
frames are written as parquet instead.

This lives in the core rather than in one subpackage: persisting a frame is
not specific to any pipeline, and a subpackage should not have to import
another one to do it.

:func:`read_frame` still accepts the old ``.p`` files, so a half-finished run
from a previous version can be picked up where it left off.
"""

import pickle
from pathlib import Path

import pandas as pd

FRAME_SUFFIX = ".parquet"
LEGACY_SUFFIX = ".p"


def _with_suffix(path: Path, suffix: str) -> Path:
    """`path` carrying `suffix`, without eating part of the name.

    ``Path.with_suffix`` replaces everything after the last dot, so a sample
    named ``pt.01_MTP`` would become ``pt.parquet`` — and two samples could
    collapse onto the same file. Only a suffix this module owns is replaced.
    """
    path = Path(path)
    if path.suffix in (FRAME_SUFFIX, LEGACY_SUFFIX):
        return path.with_suffix(suffix)
    return path.with_name(path.name + suffix)


def frame_path(path: Path) -> Path:
    """The parquet path for `path`, whatever suffix it was given."""
    return _with_suffix(path, FRAME_SUFFIX)


def legacy_path(path: Path) -> Path:
    """The pre-parquet pickle path for `path`."""
    return _with_suffix(path, LEGACY_SUFFIX)


def write_frame(df: pd.DataFrame, path: Path) -> Path:
    """Write `df` as parquet, returning the path actually written."""
    target = frame_path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    # pyarrow is a core dependency and handles the list-valued columns the
    # detection and validation stages produce.
    df.to_parquet(target, engine="pyarrow")
    return target


def read_frame(path: Path) -> pd.DataFrame:
    """Read a frame written by :func:`write_frame`, or a legacy pickle."""
    target = frame_path(path)
    if target.exists():
        # The default engine resolves to pyarrow, which is a core dependency;
        # naming it explicitly trips a pandas-stubs overload.
        return pd.read_parquet(target)

    legacy = legacy_path(path)
    if legacy.exists():
        with open(legacy, "rb") as f:
            return pickle.load(f)

    raise FileNotFoundError(f"No frame at {target} (or {legacy})")


def frame_exists(path: Path) -> bool:
    """Whether a frame is available at `path`, in either format."""
    return frame_path(path).exists() or legacy_path(path).exists()
