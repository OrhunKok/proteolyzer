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


@pytest.fixture
def jmod_ids() -> pd.DataFrame:
    """A JMod identification table, under JMod's own column names."""
    rows = 6
    return pd.DataFrame(
        {
            "file_name": ["run1", "run2"] * (rows // 2),
            "seq": [f"PEPTIDEK{i}(mTRAQ-K-4)" for i in range(rows)],
            "stripped_seq": [f"PEPTIDEK{i}" for i in range(rows)],
            "z": [2, 3] * (rows // 2),
            "silac_channel": ["mTRAQ-K-4"] * rows,
            "protein": [f"P{i:05d}" for i in range(rows)],
            "rt": np.linspace(10.0, 40.0, rows),
            "mz": np.linspace(400.0, 900.0, rows),
            "plex_Area": np.linspace(1e5, 1e7, rows),
            "MS1_Area": np.linspace(1e5, 1e7, rows),
            "Qvalue": np.linspace(0.001, 0.01, rows),
            "pep_len": [8] * rows,
            # Something for a cols_to_load test to leave behind.
            "unused_column": np.arange(rows, dtype=float),
        }
    )


@pytest.fixture
def spectronaut_report() -> pd.DataFrame:
    """A Spectronaut report, under Spectronaut's own column names.

    Long format, one row a precursor a run, with the things a real export has
    that no other format here does: names carrying a space and brackets, booleans
    where DIA-NN writes 1/0, a quantity column spanning three orders of magnitude
    from 2.54 up, and the isolation window each precursor was taken in. No
    ``EG.PrecursorId``, as the export this was written from had none.
    """
    rows = 6
    return pd.DataFrame(
        {
            "R.FileName": ["run1", "run2"] * (rows // 2),
            "PG.ProteinGroups": [f"P{i:05d}" for i in range(rows)],
            "PG.Quantity": np.linspace(1e5, 1e7, rows),
            "PG.Qvalue": np.linspace(0.0001, 0.001, rows),
            "PG.Cscore (Run-Wise)": np.linspace(1.0, 4.0, rows),
            "PEP.StrippedSequence": [f"PEPTIDEK{i}" for i in range(rows)],
            "PEP.NrOfMissedCleavages": [0, 1] * (rows // 2),
            "PEP.IsProteotypic": [True, False] * (rows // 2),
            "EG.ModifiedSequence": [f"_PEPTIDEK{i}_" for i in range(rows)],
            "EG.ApexRT": np.linspace(10.0, 40.0, rows),
            "EG.RTPredicted": np.linspace(10.5, 40.5, rows),
            "EG.Qvalue": np.linspace(0.001, 0.01, rows),
            "EG.PEP": np.linspace(0.001, 0.05, rows),
            "EG.IsDecoy": [False] * rows,
            "FG.Charge": [2, 3] * (rows // 2),
            "FG.PrecMz": np.linspace(400.0, 900.0, rows),
            # 2.54 and 400,000 in the one column, as a real export has.
            "FG.Quantity": [2.54, 43.7, 812.25, 9134.5, 87021.0, 400000.0],
            # Which of the method's isolation windows took the precursor.
            "FG.PrecWindowNumber": [1, 17, 42, 3, 58, 60],
            # Something for a cols_to_load test to leave behind.
            "EG.TotalQuantity (Settings)": np.linspace(1e4, 1e6, rows),
        }
    )


@pytest.fixture
def fragpipe_psms() -> pd.DataFrame:
    """A FragPipe/Philosopher psm.tsv, under FragPipe's own column names."""
    rows = 6
    return pd.DataFrame(
        {
            "Spectrum": [f"run1.{i:05d}.{i:05d}.2" for i in range(rows)],
            "Spectrum File": ["interact-run1.pep.xml"] * rows,
            "Peptide": [f"PEPTIDEK{i}" for i in range(rows)],
            "Modified Peptide": [f"PEPTIDEK{i}" for i in range(rows)],
            "Peptide Length": [8] * rows,
            "Charge": [2, 3] * (rows // 2),
            "Retention": np.linspace(600.0, 2400.0, rows),
            "Observed M/Z": np.linspace(400.0, 900.0, rows),
            "Hyperscore": np.linspace(10.0, 40.0, rows),
            "Intensity": np.linspace(1e5, 1e7, rows),
            "Number of Missed Cleavages": [0, 1] * (rows // 2),
            "Protein": [f"sp|P{i:05d}|GENE{i}_HUMAN" for i in range(rows)],
            "Gene": [f"GENE{i}" for i in range(rows)],
            "Unused Column": np.arange(rows, dtype=float),
        }
    )
