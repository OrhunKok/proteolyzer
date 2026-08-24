import numpy as np
import pandas as pd
import pytest

from proteolyzer.core.operations import cv, jaccard_index, per_distinct


def test_cv_matches_manual_calculation():
    data = [10.0, 12.0, 14.0]
    expected = np.std(data, ddof=1) / np.mean(data)
    assert cv(data) == pytest.approx(expected)


def test_cv_needs_a_minimum_number_of_points():
    assert np.isnan(cv([1.0, 2.0]))
    assert not np.isnan(cv([1.0, 2.0], min_datapoints=2))


def test_cv_of_zero_mean_is_nan():
    assert np.isnan(cv([-1.0, 0.0, 1.0]))


def test_cv_of_constant_data_is_zero():
    assert cv([5.0, 5.0, 5.0]) == pytest.approx(0.0)


# --------------------------------------------------------------- per_distinct

#: Element-wise transformations of a string column, of the kinds proteolyzer
#: applies to identifiers, sequences and protein groups.
TRANSFORMATIONS = {
    "split": lambda x: x.str.split(";").str[0],
    "length": lambda x: x.str.len(),
    "last character": lambda x: x.str[-1],
    "comparison": lambda x: x.str[-1] == "A",
    "count": lambda x: x.str.count("A"),
}

COLUMNS = {
    "repeated values": pd.Series(["B;A", "B;A", "C", "C"]),
    "all distinct": pd.Series(["B;A", "C"]),
    "with gaps": pd.Series(["B;A", None, "C", np.nan]),
    "empty string": pd.Series(["", "X;Y", ""]),
    "categorical": pd.Series(["B;A", "C", "B;A"], dtype="category"),
    "categorical with gaps": pd.Series(["B;A", None, "C"], dtype="category"),
    "no rows": pd.Series([], dtype=object),
    "unsorted index": pd.Series(["B;A", "C"], index=[7, 3]),
    "repeated index": pd.Series(["B;A", "C"], index=[1, 1]),
}


@pytest.mark.parametrize("column", COLUMNS.values(), ids=list(COLUMNS))
@pytest.mark.parametrize(
    "transformation", TRANSFORMATIONS.values(), ids=list(TRANSFORMATIONS)
)
def test_per_distinct_matches_applying_the_function_per_row(column, transformation):
    """The whole point: an accelerator that cannot change the answer."""
    expected = transformation(column)
    result = per_distinct(transformation)(column)

    pd.testing.assert_series_equal(result, expected)


@pytest.mark.parametrize(
    "transformation", TRANSFORMATIONS.values(), ids=list(TRANSFORMATIONS)
)
def test_per_distinct_of_an_all_empty_column_differs_only_in_the_sentinel(
    transformation,
):
    """The one documented difference: None in an all-None column becomes NaN.

    ``Series.equals`` treats the two as the same missing value, so this is as
    far as the guarantee goes; the stricter assert_series_equal would not.
    """
    column = pd.Series([None, None], dtype=object)

    result = per_distinct(transformation)(column)
    expected = transformation(column)

    assert result.equals(expected)
    assert result.dtype == expected.dtype


def test_per_distinct_calls_the_function_once_per_distinct_value():
    seen = []

    def record(values):
        seen.extend(values.tolist())
        return values.str.len()

    column = pd.Series(["AAA", "AAA", "AAA", "BB"])
    per_distinct(record)(column)

    assert seen == ["AAA", "BB"]


def test_jaccard_index_is_the_shared_over_the_union():
    left = [True, True, True, False]
    right = [True, True, False, False]

    assert jaccard_index(left, right) == pytest.approx(2 / 3)


def test_jaccard_index_of_identical_masks_is_one():
    mask = [True, False, True]

    assert jaccard_index(mask, mask) == 1.0


def test_jaccard_index_of_disjoint_masks_is_zero():
    assert jaccard_index([True, False], [False, True]) == 0.0


def test_jaccard_index_of_two_empty_masks_is_nan():
    """0/0: there is nothing for them to be similar about, and nan aggregates."""
    assert np.isnan(jaccard_index([False, False], [False, False]))


def test_jaccard_index_reads_counts_as_presence():
    """The masks come off columns that may be counts rather than flags, and the
    question asked of them is which precursors are there at all."""
    assert jaccard_index([3, 0, 7], [1, 0, 0]) == pytest.approx(0.5)
