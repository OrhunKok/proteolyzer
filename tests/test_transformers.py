import logging

import numpy as np
import pandas as pd
import pytest

from proteolyzer.core.matrix import MatrixBuilder


@pytest.fixture
def long_data() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "Precursor.Id": ["p1", "p1", "p2", "p2"],
            "Run": ["r1", "r2", "r1", "r2"],
            "Ms1.Area": [10.0, 20.0, 30.0, 0.0],
        }
    )


def test_matrix_generation_pivots(long_data):
    matrix = (
        MatrixBuilder(long_data)
        .matrix_generation("Ms1.Area", ["Precursor.Id"], ["Run"])
        .matrix
    )
    assert matrix.shape == (2, 2)
    assert matrix.loc["p1", "r2"] == 20.0


def test_duplicate_combinations_are_refused(long_data):
    duplicated = pd.concat([long_data, long_data.head(1)], ignore_index=True)
    builder = MatrixBuilder(duplicated)
    with pytest.raises(ValueError, match="Duplicate combinations"):
        builder.matrix_generation("Ms1.Area", ["Precursor.Id"], ["Run"])


def test_missingness_is_reported(long_data, caplog):
    caplog.set_level(logging.INFO)
    matrix = pd.DataFrame({"r1": [1.0, np.nan], "r2": [0.0, np.nan]})
    MatrixBuilder(long_data).missingness_check(matrix)
    assert "Missing At Random" in caplog.text
    assert "Missing Not At Random" in caplog.text
    assert "recommend dropping" in caplog.text


def test_missingness_rejects_out_of_range_threshold(long_data, caplog):
    MatrixBuilder(long_data).missingness_check(
        pd.DataFrame({"r1": [1.0]}), warning_threshold=7
    )
    assert "must be between 0 and 1" in caplog.text


def test_normalize_matrix_divides_by_group_aggregate(long_data):
    builder = MatrixBuilder(long_data).matrix_generation(
        "Ms1.Area", ["Precursor.Id"], ["Run"]
    )
    builder.normalize_matrix(within_groups=["Run"], agg_func=np.nansum)
    # Each group holds a single run, so every present value normalizes to 1.
    assert builder.matrix.loc["p1", "r1"] == pytest.approx(1.0)
    # Zeros are treated as missing.
    assert np.isnan(builder.matrix.loc["p2", "r2"])


def _grouped_matrix(values) -> pd.DataFrame:
    """A matrix whose two columns belong to the same normalization group."""
    return pd.DataFrame(
        values,
        index=["p1"],
        columns=pd.MultiIndex.from_tuples(
            [("plex1", "r1"), ("plex1", "r2")], names=["Plex", "Run"]
        ),
    )


def test_normalize_matrix_aggregates_across_the_group(long_data):
    builder = MatrixBuilder(long_data)
    builder.matrix = _grouped_matrix([[1.0, 3.0]])
    builder.normalize_matrix(within_groups=["Plex"], agg_func=np.nansum)
    assert builder.matrix.iloc[0].tolist() == [0.25, 0.75]


@pytest.mark.parametrize(("replace_zeros", "expected"), [(True, np.nan), (False, 0.0)])
def test_normalize_matrix_zero_handling(long_data, replace_zeros, expected):
    builder = MatrixBuilder(long_data)
    builder.matrix = _grouped_matrix([[0.0, 4.0]])
    builder.normalize_matrix(
        within_groups=["Plex"], agg_func=np.nansum, replace_zeros=replace_zeros
    )
    first, second = builder.matrix.iloc[0].tolist()
    assert np.isnan(first) if np.isnan(expected) else first == expected
    assert second == pytest.approx(1.0)


def test_normalize_matrix_handles_nullable_integer_matrices(long_data):
    """Regression: nullable dtypes made numpy raise on pd.NA."""
    integer_data = long_data.assign(
        **{"Ms1.Area": long_data["Ms1.Area"].astype("Int64")}
    )
    builder = MatrixBuilder(integer_data).matrix_generation(
        "Ms1.Area", ["Precursor.Id"], ["Run"]
    )
    builder.normalize_matrix(within_groups=["Run"], agg_func=np.nansum)
    assert builder.matrix.to_numpy().dtype == np.float64


def test_normalize_before_generation_is_refused(long_data):
    with pytest.raises(ValueError, match="Call matrix_generation"):
        MatrixBuilder(long_data).normalize_matrix(
            within_groups=["Run"], agg_func=np.nansum
        )


def test_a_report_can_be_given_directly(report_parquet):
    """MatrixBuilder takes a Report or a plain frame."""
    import proteolyzer as pz

    report = pz.read(report_parquet).process()
    builder = MatrixBuilder(report)
    assert builder.data is report.frame

    matrix = report.matrix("Ms1.Area", ["Precursor.Id"], ["Run"]).matrix
    assert matrix.shape == (6, 2)
