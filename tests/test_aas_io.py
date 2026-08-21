"""Frames handed between stages must survive a round trip, old files included."""

import pickle

import numpy as np
import pandas as pd
import pytest

from proteolyzer.aas.io import (
    frame_exists,
    frame_path,
    read_frame,
    write_frame,
)


@pytest.fixture
def stage_frame() -> pd.DataFrame:
    """Covers every column flavour the pipeline puts in a stage output."""
    return pd.DataFrame(
        {
            "DP Base Sequence": ["AAAK", "GGGR"],
            "MS/MS scan number": [[1, 2], [3]],
            "b_ions_aas": [np.array([3, 4]), np.array([2])],
            "aa subs positional probability": [0.95, None],
            "3-frame genome substring": [False, True],
            "Charge": pd.array([2, None], dtype="Int8"),
            "Run": pd.Categorical(["r1", "r2"]),
        }
    )


def test_frames_round_trip(tmp_path, stage_frame):
    write_frame(stage_frame, tmp_path / "sample_MTP")
    restored = read_frame(tmp_path / "sample_MTP")

    assert list(restored.columns) == list(stage_frame.columns)
    assert dict(restored.dtypes) == dict(stage_frame.dtypes)
    assert restored["DP Base Sequence"].tolist() == ["AAAK", "GGGR"]
    # Sequence-valued columns come back as arrays; the values are what matter,
    # and `.isin` (how validation consumes them) accepts either.
    assert list(restored.loc[0, "MS/MS scan number"]) == [1, 2]
    assert list(restored.loc[0, "b_ions_aas"]) == [3, 4]


def test_written_path_is_parquet_whatever_suffix_is_given(tmp_path, stage_frame):
    written = write_frame(stage_frame, tmp_path / "sample_MTP.p")
    assert written == tmp_path / "sample_MTP.parquet"
    assert written.exists()


def test_parent_directories_are_created(tmp_path, stage_frame):
    written = write_frame(stage_frame, tmp_path / "MTP" / "sample_MTP")
    assert written.exists()


def test_empty_frames_round_trip(tmp_path, stage_frame):
    write_frame(stage_frame.iloc[:0], tmp_path / "empty")
    assert read_frame(tmp_path / "empty").empty


def test_legacy_pickles_are_still_readable(tmp_path, stage_frame):
    """A run half-finished under the old format can be picked up."""
    legacy = tmp_path / "sample_MTP.p"
    with open(legacy, "wb") as f:
        pickle.dump(stage_frame, f)

    assert frame_exists(tmp_path / "sample_MTP")
    pd.testing.assert_frame_equal(read_frame(tmp_path / "sample_MTP"), stage_frame)


def test_parquet_wins_over_a_stale_pickle(tmp_path, stage_frame):
    with open(tmp_path / "sample_MTP.p", "wb") as f:
        pickle.dump(stage_frame.iloc[:0], f)
    write_frame(stage_frame, tmp_path / "sample_MTP")

    assert len(read_frame(tmp_path / "sample_MTP")) == 2


def test_missing_frame_names_both_candidates(tmp_path):
    assert not frame_exists(tmp_path / "absent")
    with pytest.raises(FileNotFoundError, match="absent.parquet"):
        read_frame(tmp_path / "absent")


def test_frame_path_is_idempotent(tmp_path):
    once = frame_path(tmp_path / "x")
    assert frame_path(once) == once


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        ("sample_MTP", "sample_MTP.parquet"),
        ("sample_MTP.p", "sample_MTP.parquet"),
        ("sample_MTP.parquet", "sample_MTP.parquet"),
        # Regression: with_suffix() ate everything after the last dot, so a
        # dotted sample id lost part of its name.
        ("pt.01_MTP", "pt.01_MTP.parquet"),
        ("run_2.5_MTP", "run_2.5_MTP.parquet"),
    ],
)
def test_frame_path_only_replaces_its_own_suffixes(tmp_path, name, expected):
    assert frame_path(tmp_path / name).name == expected


def test_dotted_sample_ids_do_not_collide(tmp_path, stage_frame):
    """Regression: pt.01_MTP and pt.02_MTP both became pt.parquet."""
    first = write_frame(stage_frame, tmp_path / "pt.01_MTP")
    second = write_frame(stage_frame.iloc[:1], tmp_path / "pt.02_MTP")

    assert first != second
    assert len(read_frame(tmp_path / "pt.01_MTP")) == 2
    assert len(read_frame(tmp_path / "pt.02_MTP")) == 1


def test_dotted_legacy_pickles_are_found(tmp_path, stage_frame):
    with open(tmp_path / "pt.01_MTP.p", "wb") as f:
        pickle.dump(stage_frame, f)

    assert frame_exists(tmp_path / "pt.01_MTP")
    assert len(read_frame(tmp_path / "pt.01_MTP")) == 2
