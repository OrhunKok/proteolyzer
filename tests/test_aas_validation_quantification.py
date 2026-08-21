"""End-to-end tests for the last two stages.

Detection leaves a SAAP frame behind; validation checks it against fragment
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

BASE, SAAP_SEQ = "AAAK", "AGAK"


@pytest.fixture
def validation_inputs(aas_params):
    """The evidence/msms parquet pair and the stage-1 SAAP frame."""
    val_dir = Path(aas_params["Utils"]["Data Folder"]) / "sample_a_val"

    pd.DataFrame(
        {
            "Raw file": ["r1", "r1"],
            "Charge": [2, 2],
            "Sequence": [BASE, SAAP_SEQ],
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
                "SAAP sequence": [SAAP_SEQ],
                "SAAP position": [1],
                "aa subs": ["A to G"],
                "MS/MS scan number": [[1, 2]],
            }
        ),
        Path(aas_params["Utils"]["Output Folder"])
        / "SAAP"
        / "sample_a_SAAP_Filtered_Stage_1",
    )
    return val_dir


def test_validation_writes_both_stage_two_frames(aas_params, validation_inputs):
    Validation(aas_params).run()

    saap_dir = Path(aas_params["Utils"]["Output Folder"]) / "SAAP"
    assert (saap_dir / "sample_a_SAAP_Filtered_Stage_2.parquet").exists()
    assert (saap_dir / "sample_a_Val_Evidence_Filtered_Stage_2.parquet").exists()


def test_validation_keeps_candidates_with_fragment_support(
    aas_params, validation_inputs
):
    Validation(aas_params).run()

    validated = read_frame(
        Path(aas_params["Utils"]["Output Folder"])
        / "SAAP"
        / "sample_a_SAAP_Filtered_Stage_2"
    )
    assert validated["SAAP sequence"].tolist() == [SAAP_SEQ]
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
        / "SAAP"
        / "sample_a_SAAP_Filtered_Stage_2"
    )
    assert validated.empty


def test_validation_normalizes_label_free_intensities(aas_params, validation_inputs):
    Validation(aas_params).run()

    evidence = read_frame(
        Path(aas_params["Utils"]["Output Folder"])
        / "SAAP"
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
        n_saap = int(rng.integers(1, 12))
        sequences = [
            "".join(rng.choice(list("AGVLK"), size=int(rng.integers(3, 9))))
            for _ in range(n_saap)
        ]
        saap = pd.DataFrame(
            {
                "Raw file": rng.choice(["r1", "r2"], n_saap),
                "Charge": rng.integers(1, 4, n_saap),
                "SAAP sequence": sequences,
                "SAAP position": [int(rng.integers(0, len(s))) for s in sequences],
                "MS/MS scan number": [
                    list(
                        rng.choice(
                            np.arange(8),
                            size=int(rng.integers(1, 4)),
                            replace=False,
                        )
                    )
                    for _ in range(n_saap)
                ],
            },
            index=rng.permutation(n_saap) if trial % 3 else None,
        )
        saap["b_ions_aas"], saap["y_ions_aas"] = zip(
            *saap.apply(validation.frags_containing_aas, axis=1), strict=True
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

        row_wise = saap.apply(
            lambda row, fragments=observed: validation.frag_count(row, fragments),
            axis=1,
        ).astype("int64")

        pd.testing.assert_series_equal(
            validation.fragment_evidence(saap, observed), row_wise
        )


def test_fragment_evidence_without_any_observed_fragments(aas_params):
    validation = Validation(aas_params)
    saap = pd.DataFrame(
        {
            "Raw file": ["r1"],
            "Charge": [2],
            "SAAP sequence": [SAAP_SEQ],
            "SAAP position": [1],
            "MS/MS scan number": [[1]],
            "b_ions_aas": [np.array([2, 3])],
            "y_ions_aas": [np.array([3])],
        }
    )
    empty = pd.DataFrame(
        {key: [] for key in Validation._FRAGMENT_KEYS},
    )

    assert validation.fragment_evidence(saap, empty).tolist() == [0]


@pytest.fixture
def quantification_inputs(aas_params, validation_inputs):
    """Run validation so quantification has its stage-2 inputs."""
    Validation(aas_params).run()
    return Path(aas_params["Utils"]["Output Folder"]) / "SAAP"


def test_quantification_writes_log2_ratios(aas_params, quantification_inputs):
    Quantification(aas_params).run()

    quant = read_frame(quantification_inputs / "sample_a_SAAP_Quant")
    assert quant["SAAP sequence"].tolist() == [SAAP_SEQ]
    # Label-free: log2(250 / 1000) = -2.
    assert quant["Ratio"].tolist() == [pytest.approx(-2.0)]
    assert quant["SAAP.Sum"].tolist() == [pytest.approx(250.0)]
    assert quant["BASE.Sum"].tolist() == [pytest.approx(1000.0)]


def test_minimum_quantity_filters_low_signal_ratios(
    aas_params, quantification_inputs, capsys
):
    """Regression: the parameter was documented but never read."""
    aas_params["Quantification"]["MaxQuant"]["Minimum Quantity"] = 500

    Quantification(aas_params).run()

    quant = read_frame(quantification_inputs / "sample_a_SAAP_Quant")
    # The mistranslated peptide sums to 250, below the threshold.
    assert quant.empty


def test_minimum_quantity_of_zero_keeps_everything(aas_params, quantification_inputs):
    aas_params["Quantification"]["MaxQuant"]["Minimum Quantity"] = 0

    Quantification(aas_params).run()

    assert not read_frame(quantification_inputs / "sample_a_SAAP_Quant").empty


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
                    "SAAP sequence": [SAAP_SEQ, SAAP_SEQ],
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


def _tmt_frames():
    """A SAAP/BASE pair with two reporter channels each.

    The normalised columns carry values that would wreck the channel
    proportions if the reporter regex let them into the row sums, which is
    what its negative lookahead is for.
    """
    saap = pd.DataFrame(
        {
            "Intensity": [300.0],
            "Reporter intensity corrected 1": [30.0],
            "Reporter intensity corrected 2": [70.0],
            "Normalised Reporter intensity corrected 1": [9000.0],
            "Normalised Reporter intensity corrected 2": [11000.0],
        }
    )
    base = pd.DataFrame(
        {
            "Intensity": [1200.0],
            "Reporter intensity corrected 1": [25.0],
            "Reporter intensity corrected 2": [75.0],
            "Normalised Reporter intensity corrected 1": [4000.0],
            "Normalised Reporter intensity corrected 2": [6000.0],
        }
    )
    saap_ids = pd.DataFrame(
        {
            "DP Base Sequence": [BASE],
            "SAAP sequence": [SAAP_SEQ],
            "aa subs": ["A to G"],
        }
    )
    return saap_ids, saap, base


def test_tmt_quantification_distributes_intensity_across_channels(aas_params):
    """Each channel gets the share of the total its reporter ion accounts for."""
    saap_ids, saap, base = _tmt_frames()
    metadata = pd.DataFrame({"MQ": [1.0, 2.0]})

    quant = Quantification(aas_params)._raas(saap_ids, saap, base, metadata, "TMT")

    # 30 of 100 reporter counts, so 30% of the 300 total.
    assert quant["SAAP.Plex.1"].tolist() == [pytest.approx(90.0)]
    assert quant["SAAP.Plex.2"].tolist() == [pytest.approx(210.0)]
    assert quant["BASE.Plex.1"].tolist() == [pytest.approx(300.0)]
    assert quant["BASE.Plex.2"].tolist() == [pytest.approx(900.0)]
    # The distributed intensities sum back to the total they came from.
    assert quant["SAAP.Plex.1"][0] + quant["SAAP.Plex.2"][0] == pytest.approx(300.0)


def test_tmt_ratios_are_log10_per_channel_and_overall(aas_params):
    saap_ids, saap, base = _tmt_frames()
    metadata = pd.DataFrame({"MQ": [1.0, 2.0]})

    quant = Quantification(aas_params)._raas(saap_ids, saap, base, metadata, "TMT")

    assert quant["Ratio"].tolist() == [pytest.approx(np.log10(300 / 1200))]
    assert quant["Ratio.Plex.1"].tolist() == [pytest.approx(np.log10(90 / 300))]
    assert quant["Ratio.Plex.2"].tolist() == [pytest.approx(np.log10(210 / 900))]


def test_tmt_normalised_reporters_are_passed_through(aas_params):
    """They are already normalised, so they are reported, not redistributed."""
    saap_ids, saap, base = _tmt_frames()
    metadata = pd.DataFrame({"MQ": [1.0, 2.0]})

    quant = Quantification(aas_params)._raas(saap_ids, saap, base, metadata, "TMT")

    assert quant["SAAP.Plex.1.Norm.Sum"].tolist() == [pytest.approx(9000.0)]
    assert quant["BASE.Plex.2.Norm.Sum"].tolist() == [pytest.approx(6000.0)]


def test_tmt_only_reports_the_channels_the_metadata_names(aas_params):
    """A plex the sample was not labelled with gets no columns."""
    saap_ids, saap, base = _tmt_frames()
    metadata = pd.DataFrame({"MQ": [2.0, np.nan]})

    quant = Quantification(aas_params)._raas(saap_ids, saap, base, metadata, "TMT")

    assert "SAAP.Plex.2" in quant.columns
    assert "SAAP.Plex.1" not in quant.columns


def test_an_unknown_label_designation_is_refused(aas_params):
    saap_ids, saap, base = _tmt_frames()
    with pytest.raises(ValueError, match="Unknown label designation"):
        Quantification(aas_params)._raas(
            saap_ids, saap, base, pd.DataFrame({"MQ": [1.0]}), "iTRAQ"
        )


def test_results_read_back_a_real_run(aas_params, quantification_inputs):
    """The other side of the pipeline: find what the stages just wrote."""
    from proteolyzer.aas.results import Results

    Quantification(aas_params).run()
    results = Results.from_params(aas_params)

    assert results.samples == ["sample_a"]
    assert results.has("quantified", "sample_a")

    summary = results.summary()
    assert summary.loc["sample_a", "validated"] == 1
    assert summary.loc["sample_a", "quantified"] == 1
    # Detection never ran here, so its artefacts are absent rather than empty.
    assert pd.isna(summary.loc["sample_a", "candidates"])

    combined = results.combined("quantified")
    assert combined["Sample"].tolist() == ["sample_a"]
    assert combined["Ratio"].tolist() == [pytest.approx(-2.0)]

    assert results.provenance()["step"].tolist() == ["Validation", "Quantification"]
