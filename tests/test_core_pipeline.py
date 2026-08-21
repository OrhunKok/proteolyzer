"""Stage plumbing and the provenance log.

This was covered only through the pipelines that used to live in this
repository. It is the surface they still depend on from their own repos, so it
is tested here directly now.
"""

import json

import pytest

from proteolyzer.core.pipeline import (
    PROVENANCE_FILE,
    NullQueue,
    Stage,
    package_version,
    record_run,
)


class Recorder:
    """A queue that keeps what it was given."""

    def __init__(self):
        self.items = []

    def put(self, item):
        self.items.append(item)


def test_a_stage_without_a_queue_discards_its_progress():
    """A stage run outside multiprocessing should not need a queue."""
    stage = Stage({"Utils": {}})

    assert isinstance(stage.queue, NullQueue)
    stage.announce()  # would raise if NullQueue did not accept it


def test_a_stage_announces_itself_by_class_name():
    class Detection(Stage):
        pass

    queue = Recorder()
    Detection({}, queue).announce()

    assert queue.items == [("stdout", "Detection initialized.")]


def test_run_is_left_to_the_subclass():
    with pytest.raises(NotImplementedError):
        Stage({}).run()


def test_a_stage_records_itself_under_its_own_name(tmp_path):
    class Quantification(Stage):
        pass

    log = Quantification({"Minimum Quantity": 0}).record_run(tmp_path)

    entry = json.loads(log.read_text())
    assert log.name == PROVENANCE_FILE
    assert entry["step"] == "Quantification"
    assert entry["params"] == {"Minimum Quantity": 0}


def test_provenance_accumulates_one_line_per_run(tmp_path):
    """Appended, not overwritten: the log is the history of the folder."""
    record_run(tmp_path, "Detection", {"a": 1})
    record_run(tmp_path, "Validation", {"a": 2})

    lines = (tmp_path / PROVENANCE_FILE).read_text().splitlines()
    assert [json.loads(line)["step"] for line in lines] == ["Detection", "Validation"]


def test_a_record_says_what_ran_when_and_from_which_version(tmp_path):
    log = record_run(tmp_path, "Detection", {})

    entry = json.loads(log.read_text())
    assert set(entry) == {"step", "timestamp", "proteolyzer_version", "params"}
    assert entry["proteolyzer_version"] == package_version()
    # An ISO timestamp to the second, in UTC.
    assert entry["timestamp"].endswith("+00:00")


def test_extra_sections_are_merged_into_the_record(tmp_path):
    log = record_run(tmp_path, "CoordinatesMapping", {}, inputs={"Label": 4}, clashes=0)

    entry = json.loads(log.read_text())
    assert entry["inputs"] == {"Label": 4}
    assert entry["clashes"] == 0


def test_the_output_directory_is_created_if_absent(tmp_path):
    target = tmp_path / "results" / "run1"

    log = record_run(target, "Detection", {})

    assert log.exists()
    assert log.parent == target


def test_unserializable_parameters_are_rendered_rather_than_refused(tmp_path):
    """A parameter schema hands over Path values, which JSON cannot encode."""
    log = record_run(tmp_path, "Detection", {"Genome FASTA": tmp_path / "genome.fa"})

    entry = json.loads(log.read_text())
    assert entry["params"]["Genome FASTA"] == str(tmp_path / "genome.fa")


def test_the_version_is_reported_rather_than_raising():
    assert isinstance(package_version(), str)
    assert package_version()
