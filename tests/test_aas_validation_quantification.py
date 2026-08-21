"""End-to-end tests for the last two stages.

Detection leaves an MTP frame behind; validation checks it against fragment
ions and quantification turns what survives into ratios. Both are driven here
from hand-built stage inputs so the assertions stay readable.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("yaml")
pytest.importorskip("fastparquet")

from proteolyzer.aas.quantification import Quantification  # noqa: E402
from proteolyzer.aas.validation import Validation  # noqa: E402
from proteolyzer.core.io import read_frame, write_frame  # noqa: E402

BASE, MTP_SEQ = "AAAK", "AGAK"


@pytest.fixture
def validation_inputs(aas_params):
    """The evidence/msms parquet pair and the stage-1 MTP frame."""
    val_dir = Path(aas_params["Utils"]["Data Folder"]) / "sample_a_val"

    pd.DataFrame(
        {
            "Raw file": ["r1", "r1"],
            "Charge": [2, 2],
            "Sequence": [BASE, MTP_SEQ],
            "MS/MS scan number": [1, 2],
            "PEP": [0.001, 0.001],
            "PIF": [0.95, 0.95],
            "Intensity": [1000.0, 250.0],
            "m/z": [500.0, 500.0],
            "Mass": [1000.0, 1000.0],
            "Retention time": [10.0, 10.5],
            "Mass error [ppm]": [1.0, 2.0],
        }
    ).to_parquet(val_dir / "evidence.parquet", engine="fastparquet")

    # Two b ions and two y ions per scan, spanning the substituted residue.
    pd.DataFrame(
        {
            "Raw file": ["r1"] * 8,
            "MS/MS scan number": [1, 1, 1, 1, 2, 2, 2, 2],
            "Frag.Type": ["b", "b", "y", "y"] * 2,
            "Frag.Number": [2, 3, 2, 3] * 2,
        }
    ).to_parquet(val_dir / "msms.parquet", engine="fastparquet")

    write_frame(
        pd.DataFrame(
            {
                "Raw file": ["r1"],
                "Charge": [2],
                "DP Base Sequence": [BASE],
                "mistranslated sequence": [MTP_SEQ],
                "mistranslated aas positions": [1],
                "aa subs": ["A to G"],
                "MS/MS scan number": [[1, 2]],
            }
        ),
        Path(aas_params["Utils"]["Output Folder"])
        / "MTP"
        / "sample_a_MTP_Filtered_Stage_1",
    )
    return val_dir


def test_validation_writes_both_stage_two_frames(aas_params, validation_inputs):
    Validation(aas_params).run()

    mtp_dir = Path(aas_params["Utils"]["Output Folder"]) / "MTP"
    assert (mtp_dir / "sample_a_MTP_Filtered_Stage_2.parquet").exists()
    assert (mtp_dir / "sample_a_Val_Evidence_Filtered_Stage_2.parquet").exists()


def test_validation_keeps_candidates_with_fragment_support(
    aas_params, validation_inputs
):
    Validation(aas_params).run()

    validated = read_frame(
        Path(aas_params["Utils"]["Output Folder"])
        / "MTP"
        / "sample_a_MTP_Filtered_Stage_2"
    )
    assert validated["mistranslated sequence"].tolist() == [MTP_SEQ]
    # Only the scan whose evidence row is the mistranslated peptide counts
    # (scan 2). Of its four fragments, b2, b3 and y3 span position 1 of AGAK.
    assert validated["fragment_evidence"].tolist() == [3]


def test_validation_drops_candidates_without_fragment_support(
    aas_params, validation_inputs
):
    """Fragments that do not span the substitution are no evidence for it."""
    pd.DataFrame(
        {
            "Raw file": ["r1"],
            "MS/MS scan number": [1],
            "Frag.Type": ["b"],
            "Frag.Number": [1],
        }
    ).to_parquet(validation_inputs / "msms.parquet", engine="fastparquet")

    Validation(aas_params).run()

    validated = read_frame(
        Path(aas_params["Utils"]["Output Folder"])
        / "MTP"
        / "sample_a_MTP_Filtered_Stage_2"
    )
    assert validated.empty


def test_validation_normalizes_label_free_intensities(aas_params, validation_inputs):
    Validation(aas_params).run()

    evidence = read_frame(
        Path(aas_params["Utils"]["Output Folder"])
        / "MTP"
        / "sample_a_Val_Evidence_Filtered_Stage_2"
    )
    assert "Intensity.Normalised" in evidence.columns


def test_fragment_evidence_matches_the_row_wise_count(aas_params):
    """The vectorized count must agree with `frag_count` exactly.

    Randomized rather than hand-picked, because the two disagree only in
    corners: shuffled indices, repeated sequences, no observed fragments, and
    ion types (a ions) that must not be counted.
    """
    validation = Validation(aas_params)
    rng = np.random.default_rng(0)

    for trial in range(25):
        n_mtp = int(rng.integers(1, 12))
        sequences = [
            "".join(rng.choice(list("AGVLK"), size=int(rng.integers(3, 9))))
            for _ in range(n_mtp)
        ]
        mtp = pd.DataFrame(
            {
                "Raw file": rng.choice(["r1", "r2"], n_mtp),
                "Charge": rng.integers(1, 4, n_mtp),
                "mistranslated sequence": sequences,
                "mistranslated aas positions": [
                    int(rng.integers(0, len(s))) for s in sequences
                ],
                "MS/MS scan number": [
                    list(
                        rng.choice(
                            np.arange(8),
                            size=int(rng.integers(1, 4)),
                            replace=False,
                        )
                    )
                    for _ in range(n_mtp)
                ],
            },
            index=rng.permutation(n_mtp) if trial % 3 else None,
        )
        mtp["b_ions_aas"], mtp["y_ions_aas"] = zip(
            *mtp.apply(validation.frags_containing_aas, axis=1), strict=True
        )

        n_msms = int(rng.integers(0, 60))
        observed = pd.DataFrame(
            {
                "Raw file": rng.choice(["r1", "r2"], n_msms),
                "Sequence": rng.choice([*sequences, "ZZZZ"], n_msms),
                "MS/MS scan number": rng.integers(0, 8, n_msms),
                "Charge": rng.integers(1, 4, n_msms),
                "Frag.Type": rng.choice(["a", "b", "y"], n_msms),
                "Frag.Number": rng.integers(0, 9, n_msms),
            }
        )

        row_wise = mtp.apply(
            lambda row, fragments=observed: validation.frag_count(row, fragments),
            axis=1,
        ).astype("int64")

        pd.testing.assert_series_equal(
            validation.fragment_evidence(mtp, observed), row_wise
        )


def test_fragment_evidence_without_any_observed_fragments(aas_params):
    validation = Validation(aas_params)
    mtp = pd.DataFrame(
        {
            "Raw file": ["r1"],
            "Charge": [2],
            "mistranslated sequence": [MTP_SEQ],
            "mistranslated aas positions": [1],
            "MS/MS scan number": [[1]],
            "b_ions_aas": [np.array([2, 3])],
            "y_ions_aas": [np.array([3])],
        }
    )
    empty = pd.DataFrame(
        {key: [] for key in Validation._FRAGMENT_KEYS},
    )

    assert validation.fragment_evidence(mtp, empty).tolist() == [0]


@pytest.fixture
def quantification_inputs(aas_params, validation_inputs):
    """Run validation so quantification has its stage-2 inputs."""
    Validation(aas_params).run()
    return Path(aas_params["Utils"]["Output Folder"]) / "MTP"


def test_quantification_writes_log2_ratios(aas_params, quantification_inputs):
    Quantification(aas_params).run()

    quant = read_frame(quantification_inputs / "sample_a_MTP_Quant")
    assert quant["mistranslated sequence"].tolist() == [MTP_SEQ]
    # Label-free: log2(250 / 1000) = -2.
    assert quant["Ratio"].tolist() == [pytest.approx(-2.0)]
    assert quant["MTP.Sum"].tolist() == [pytest.approx(250.0)]
    assert quant["BP.Sum"].tolist() == [pytest.approx(1000.0)]


def test_minimum_quantity_filters_low_signal_ratios(
    aas_params, quantification_inputs, capsys
):
    """Regression: the parameter was documented but never read."""
    aas_params["Quantification"]["MaxQuant"]["Minimum Quantity"] = 500

    Quantification(aas_params).run()

    quant = read_frame(quantification_inputs / "sample_a_MTP_Quant")
    # The mistranslated peptide sums to 250, below the threshold.
    assert quant.empty


def test_minimum_quantity_of_zero_keeps_everything(aas_params, quantification_inputs):
    aas_params["Quantification"]["MaxQuant"]["Minimum Quantity"] = 0

    Quantification(aas_params).run()

    assert not read_frame(quantification_inputs / "sample_a_MTP_Quant").empty


def test_quantification_reports_missing_validation_directory(aas_params):
    class RecordingQueue:
        def __init__(self):
            self.messages = []

        def put(self, item):
            self.messages.append(item)

    queue = RecordingQueue()
    Quantification(aas_params, queue).process_sample("sample_b")
    assert ("stderr", "sample_b validation directory not found") in queue.messages


def test_quantification_ratio_is_finite_only(aas_params, quantification_inputs):
    """A zero base-peptide intensity yields -inf, which must not be written."""
    with np.errstate(divide="ignore"):
        quant = Quantification(aas_params)._raas(
            pd.DataFrame(
                {
                    "DP Base Sequence": [BASE, BASE],
                    "mistranslated sequence": [MTP_SEQ, MTP_SEQ],
                    "aa subs": ["A to G", "A to G"],
                }
            ),
            pd.DataFrame({"Intensity": [250.0, 250.0]}),
            pd.DataFrame({"Intensity": [1000.0, 0.0]}),
            pd.DataFrame({"MQ": [np.nan]}),
            "Label-Free",
        )
    assert np.isinf(quant["Ratio"]).any()
    assert len(quant[np.isfinite(quant["Ratio"])]) == 1
