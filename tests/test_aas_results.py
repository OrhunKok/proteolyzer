"""Reading back an AAS output folder.

The stage tests check that the frames are written; these check that someone
can find them afterwards without knowing the stage-internal file names.
"""

from pathlib import Path

import pandas as pd
import pytest

pytest.importorskip("yaml")

from proteolyzer.aas.results import ARTEFACTS, SAMPLE_COL, Results  # noqa: E402
from proteolyzer.core.io import write_frame  # noqa: E402
from proteolyzer.core.pipeline import record_run  # noqa: E402


def _write(results: Results, artefact: str, sample: str, rows: int) -> None:
    """Put `rows` rows of `artefact` where the pipeline would have."""
    frame = pd.DataFrame(
        {
            "DP Base Sequence": ["AAAK"] * rows,
            "aa subs": ["A to G"] * rows,
            "Ratio": [-2.0] * rows,
        }
    )
    write_frame(frame, results.path(artefact, sample))


@pytest.fixture
def results(tmp_path) -> Results:
    """A folder where two samples got as far as different steps."""
    results = Results(tmp_path / "out")
    _write(results, "candidates", "sample_a", 5)
    _write(results, "validated", "sample_a", 3)
    _write(results, "quantified", "sample_a", 2)
    # sample_b stopped after detection.
    _write(results, "candidates", "sample_b", 4)
    return results


def test_samples_are_discovered(results):
    assert results.samples == ["sample_a", "sample_b"]


def test_sample_ids_containing_underscores_survive(tmp_path):
    """Splitting on '_' rather than stripping the stem would truncate these."""
    results = Results(tmp_path / "out")
    _write(results, "quantified", "kelly_single_cell_012126", 1)
    assert results.samples == ["kelly_single_cell_012126"]


def test_an_empty_folder_has_no_samples(tmp_path):
    assert Results(tmp_path / "nothing").samples == []


def test_which_steps_reached_which_sample(results):
    assert results.has("quantified", "sample_a")
    assert not results.has("quantified", "sample_b")


def test_an_unknown_artefact_lists_the_real_ones(results):
    with pytest.raises(ValueError, match="No artefact 'nope'"):
        results.path("nope", "sample_a")


def test_one_sample_is_loaded_by_result_name(results):
    """The caller names the result, not the stage's file."""
    assert len(results.load("quantified", "sample_a")) == 2


def test_samples_are_combined_with_a_sample_column(results):
    combined = results.combined("candidates")

    assert list(combined.columns)[0] == SAMPLE_COL
    assert combined[SAMPLE_COL].tolist() == ["sample_a"] * 5 + ["sample_b"] * 4
    assert len(combined) == 9


def test_combining_reports_the_samples_a_step_did_not_reach(results, caplog):
    combined = results.combined("quantified")

    assert combined[SAMPLE_COL].unique().tolist() == ["sample_a"]
    assert "No 'quantified' for 1 of 2 samples: ['sample_b']" in caplog.text


def test_combining_nothing_gives_an_empty_frame(results):
    combined = results.combined("evidence")
    assert combined.empty
    assert list(combined.columns) == [SAMPLE_COL]


def test_a_subset_of_samples_can_be_asked_for(results):
    combined = results.combined("candidates", samples=["sample_b"])
    assert combined[SAMPLE_COL].unique().tolist() == ["sample_b"]


def test_the_summary_shows_where_the_pipeline_stopped(results):
    summary = results.summary()

    assert list(summary.index) == ["sample_a", "sample_b"]
    assert list(summary.columns) == list(ARTEFACTS)
    assert summary.loc["sample_a", "quantified"] == 2
    assert pd.isna(summary.loc["sample_b", "quantified"])
    assert summary.loc["sample_b", "candidates"] == 4


def test_provenance_is_read_back(results):
    record_run(results.output_dir, "Detection", {"Utils": {"Label Plex": 1}})
    record_run(results.output_dir, "Validation", {"Utils": {"Label Plex": 1}})

    provenance = results.provenance()

    assert provenance["step"].tolist() == ["Detection", "Validation"]
    assert provenance["proteolyzer_version"].notna().all()


def test_a_folder_without_provenance_is_reported(results, caplog):
    assert results.provenance().empty
    assert "No provenance.jsonl" in caplog.text


def test_results_can_be_opened_from_the_parameters(aas_params):
    results = Results.from_params(aas_params)
    assert results.output_dir == Path(aas_params["Utils"]["Output Folder"])


def test_a_folder_from_before_the_saap_naming_is_still_read(tmp_path):
    """Runs written as MTP/ and PTM/ must not become unreadable."""
    results = Results(tmp_path / "out")
    frame = pd.DataFrame({"aa subs": ["A to G"], "Ratio": [-2.0]})
    write_frame(frame, tmp_path / "out" / "MTP" / "old_sample_MTP_Quant")
    write_frame(frame, tmp_path / "out" / "PTM" / "old_sample_PTM")

    assert results.samples == ["old_sample"]
    assert results.has("quantified", "old_sample")
    assert results.has("alt", "old_sample")
    assert len(results.load("quantified", "old_sample")) == 1
    assert results.combined("quantified")[SAMPLE_COL].tolist() == ["old_sample"]


def test_the_current_layout_wins_over_a_legacy_copy(tmp_path):
    results = Results(tmp_path / "out")
    write_frame(
        pd.DataFrame({"n": [1, 2, 3]}), tmp_path / "out" / "MTP" / "s_MTP_Quant"
    )
    write_frame(pd.DataFrame({"n": [1]}), tmp_path / "out" / "SAAP" / "s_SAAP_Quant")

    assert len(results.load("quantified", "s")) == 1
