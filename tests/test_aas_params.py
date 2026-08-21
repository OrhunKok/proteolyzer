from pathlib import Path

import pytest

yaml = pytest.importorskip("yaml")

from proteolyzer.aas.params import load_params  # noqa: E402


def test_workflow_sections_are_merged_into_their_parent(aas_params):
    params = load_params(aas_params)
    detection = params["Detection"]
    # Shared keys stay, and the MaxQuant block is flattened in beside them.
    assert detection["Protease"] == "Trypsin"
    assert detection["PTM ppm"] == 5.0
    assert "MaxQuant" not in detection


def test_paths_are_parsed_as_paths(aas_params):
    params = load_params(aas_params)
    assert isinstance(params["Utils"]["Output Folder"], Path)
    assert isinstance(params["Translation"]["Genome FASTA"], Path)


def test_params_can_be_loaded_from_yaml(tmp_path, aas_params):
    path = tmp_path / "params.yaml"
    path.write_text(yaml.safe_dump(aas_params))
    assert load_params(path) == load_params(str(path)) == load_params(aas_params)


def test_missing_file_is_reported(tmp_path):
    with pytest.raises(ValueError, match="does not exist"):
        load_params(tmp_path / "absent.yaml")


def test_wrong_argument_type_is_reported():
    with pytest.raises(TypeError, match="file path or dict"):
        load_params(3)


def test_unknown_keys_are_rejected(aas_params):
    aas_params["Utils"]["Typo"] = 1
    with pytest.raises(ValueError, match="Typo"):
        load_params(aas_params)


def test_unsupported_workflow_is_rejected(aas_params):
    aas_params["Utils"]["Workflow"] = "Spectronaut"
    with pytest.raises(ValueError, match="Workflow"):
        load_params(aas_params)


def test_missing_section_is_rejected(aas_params):
    del aas_params["Validation"]
    with pytest.raises(ValueError, match="Validation"):
        load_params(aas_params)


def test_shipped_example_params_are_valid_apart_from_paths():
    """The example file documents the schema, so it must stay in sync with it."""
    example = Path(__file__).resolve().parents[1] / "examples" / "aas" / "params.yaml"
    raw = yaml.safe_load(example.read_text())
    assert set(raw) == {
        "Utils",
        "Translation",
        "Detection",
        "Validation",
        "Quantification",
    }
    # Placeholders are paths, so the schema accepts them as-is.
    assert load_params(raw)["Utils"]["Workflow"] == "MaxQuant"
