"""Shared fixtures.

Everything here is synthetic: the tests must not depend on real search output,
so each fixture builds the smallest frame that still exercises the code path.
"""

import numpy as np
import pandas as pd
import pytest

#: Columns DIA-NN writes that are *not* in the configured load list. Tests use
#: them to check that column subsetting actually happens.
EXTRA_REPORT_COLS = ("Q.Value", "Predicted.RT")


def _report_frame(precursor_ids: list[str], runs: list[str]) -> pd.DataFrame:
    """A DIA-NN-shaped precursor report."""
    n = len(precursor_ids)
    stripped = [
        pid.split("(")[0].rstrip("0123456789") or "AAAK" for pid in precursor_ids
    ]
    return pd.DataFrame(
        {
            "Run": runs,
            "Precursor.Id": precursor_ids,
            "Stripped.Sequence": stripped,
            "Precursor.Charge": [2, 3] * (n // 2) + [2] * (n % 2),
            "Protein.Group": [f"P{i:05d}" for i in range(n)],
            "Genes": [f"GENE{i}" for i in range(n)],
            "Proteotypic": [1, 0] * (n // 2) + [1] * (n % 2),
            "RT": np.linspace(10.0, 40.0, n),
            "RT.Start": np.linspace(9.5, 39.5, n),
            "RT.Stop": np.linspace(10.5, 40.5, n),
            "PEP": np.linspace(0.001, 0.05, n),
            "Ms1.Area": np.linspace(1e5, 1e7, n),
            "Precursor.Quantity": np.linspace(2e5, 2e7, n),
            **{col: np.linspace(0.0001, 0.01, n) for col in EXTRA_REPORT_COLS},
        }
    )


@pytest.fixture
def report_frame():
    """Build a report frame from explicit precursor ids and runs.

    For tests that need a shape the fixtures below do not have -- a precursor
    repeated across runs, say.
    """
    return _report_frame


@pytest.fixture
def label_free_report() -> pd.DataFrame:
    ids = ["AAAKPEPTIDER2", "MYSEQK2", "VLDATRK3", "GGGGR2", "SAMPLEK2", "TESTKR3"]
    runs = ["run1", "run2", "run1", "run2", "run1", "run2"]
    return _report_frame(ids, runs)


@pytest.fixture
def labelled_report() -> pd.DataFrame:
    """A labelled report whose last precursor carries no label.

    Precursors without a channel are what used to break run-channel generation.
    """
    ids = [
        "AAAKPEPTIDER(mTRAQ-R-0)2",
        "MYSEQK(mTRAQ-K-8)2",
        "K(mTRAQ-K-0)VLDATR(mTRAQ-R-0)3",
        "GGGGR(mTRAQ-R-8)2",
        "SAMPLEK(mTRAQ-K-8)2",
        "UNLABELLEDPEPTIDEK2",
    ]
    runs = ["run1", "run2", "run1", "run2", "run1", "run2"]
    return _report_frame(ids, runs)


@pytest.fixture
def report_parquet(tmp_path, label_free_report):
    """`label_free_report` written where DIA-NN would put it."""
    path = tmp_path / "report.parquet"
    label_free_report.to_parquet(path)
    return path
