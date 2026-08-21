"""The two-phase orchestrator.

The stages themselves are tested elsewhere; these check the ordering the
orchestrator exists for, with the stages replaced by recorders.
"""

import pickle
from pathlib import Path

import pytest

pytest.importorskip("yaml")

from proteolyzer.aas.pipeline import Pipeline  # noqa: E402


@pytest.fixture
def recorded(monkeypatch) -> list[str]:
    """Replace every stage with something that just records that it ran."""
    ran: list[str] = []

    def recorder(name):
        class Recorder:
            def __init__(self, params, queue=None):
                pass

            def run(self):
                ran.append(name)

        return Recorder

    monkeypatch.setattr(
        "proteolyzer.aas.pipeline.Preprocessor",
        type("P", (), {"MaxQuant": recorder("preprocess")}),
    )
    for attr, name in (
        ("FrameTranslator", "translate"),
        ("Detection", "detect"),
        ("Validation", "validate"),
        ("Quantification", "quantify"),
    ):
        monkeypatch.setattr(f"proteolyzer.aas.pipeline.{attr}", recorder(name))
    return ran


def _write_frames(aas_params) -> None:
    """Stand in for a translation that has already happened."""
    frames = Path(aas_params["Translation"]["Translated Frames Folder"])
    for frame in range(1, 7):
        (frames / f"frame_{frame}.p").write_bytes(pickle.dumps("MTTT*"))


def _make_validation_search(aas_params, name="sample_a_val") -> Path:
    search = Path(aas_params["Utils"]["Data Folder"]) / name
    search.mkdir(exist_ok=True)
    return search


def test_phase_one_runs_in_order(aas_params, recorded):
    Pipeline(aas_params).run_detection()
    assert recorded == ["preprocess", "translate", "detect"]


def test_translation_is_skipped_when_its_frames_exist(aas_params, recorded, caplog):
    """It is the slow step and its output only depends on the genome."""
    _write_frames(aas_params)

    Pipeline(aas_params).run_detection()

    assert recorded == ["preprocess", "detect"]
    assert "skipping translation" in caplog.text


def test_translation_can_be_forced(aas_params, recorded):
    _write_frames(aas_params)
    Pipeline(aas_params).run_detection(translate=True)
    assert recorded == ["preprocess", "translate", "detect"]


def test_translation_can_be_skipped_explicitly(aas_params, recorded):
    Pipeline(aas_params).run_detection(translate=False)
    assert recorded == ["preprocess", "detect"]


def test_phase_two_preprocesses_again_before_validating(aas_params, recorded):
    """The validation searches did not exist when phase one preprocessed."""
    Pipeline(aas_params).run_validation()
    assert recorded == ["preprocess", "validate", "quantify"]


def test_phase_two_refuses_to_run_before_the_search(aas_params, recorded, tmp_path):
    """Without the manual search there is nothing for validation to read."""
    empty = tmp_path / "no_searches"
    empty.mkdir()
    aas_params["Utils"]["Data Folder"] = str(empty)

    with pytest.raises(FileNotFoundError, match="search of the validation FASTA"):
        Pipeline(aas_params).run_validation()
    assert recorded == []


def test_status_says_what_can_run(aas_params):
    pipeline = Pipeline(aas_params)

    status = pipeline.status()
    assert status["frames_translated"] is False
    assert status["can_run_detection"] is True
    # The fixture already has a sample_a_val directory.
    assert status["validation_searches"] == ["sample_a_val"]
    assert status["can_run_validation"] is True

    _write_frames(aas_params)
    assert Pipeline(aas_params).status()["frames_translated"] is True


def test_status_reports_no_searches_yet(aas_params, tmp_path):
    empty = tmp_path / "no_searches"
    empty.mkdir()
    aas_params["Utils"]["Data Folder"] = str(empty)

    status = Pipeline(aas_params).status()
    assert status["validation_searches"] == []
    assert status["can_run_validation"] is False


def test_the_phases_return_the_results(aas_params, recorded):
    pipeline = Pipeline(aas_params)
    assert pipeline.run_detection().output_dir == pipeline.output_dir
    assert pipeline.run_validation().output_dir == pipeline.output_dir


def test_extra_validation_searches_are_found(aas_params, recorded):
    _make_validation_search(aas_params, "sample_b_val")
    status = Pipeline(aas_params).status()
    assert status["validation_searches"] == ["sample_a_val", "sample_b_val"]
