import numpy as np
import pytest

from proteolyzer.utils.operations import cv


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
