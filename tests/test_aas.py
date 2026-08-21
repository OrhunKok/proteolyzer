"""Tests for the AAS reference tables and the shared stage plumbing."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("yaml")

from proteolyzer.aas.base import PROVENANCE_FILE, NullQueue, Stage  # noqa: E402
from proteolyzer.aas.io import read_frame  # noqa: E402
from proteolyzer.aas.utils import (  # noqa: E402
    aa_subs_ref,
    calculate_aa_substitution_matrix,
    column_mapping,
    gen_mod_dict,
    ptm_mtp_output,
)


class RecordingQueue:
    def __init__(self):
        self.messages = []

    def put(self, item):
        self.messages.append(item)


# --------------------------------------------------------------------- tables


def test_aa_subs_ref_is_keyed_by_origin_residue():
    subs = aa_subs_ref()
    assert subs["A"]["A to G"] == pytest.approx(-14.01565, abs=1e-5)
    # Every entry is a substitution away from its key.
    assert all(
        name.startswith(f"{origin} to ")
        for origin, table in subs.items()
        for name in table
    )


def test_aa_subs_ref_masses_match_the_amino_acid_table():
    from proteolyzer import config

    subs = aa_subs_ref()
    expected = config.AminoAcids.MASS["G"] - config.AminoAcids.MASS["A"]
    assert subs["A"]["A to G"] == pytest.approx(expected, abs=1e-5)


def test_gen_mod_dict_rows_are_name_position_mass():
    mods = gen_mod_dict()
    name, position, mass = mods["A"][0]
    assert isinstance(name, str)
    assert isinstance(position, str)
    assert isinstance(mass, float)


def test_gen_mod_dict_excludes_amino_acid_substitutions():
    """The PTM reference must not contain the substitutions it is compared against.

    Otherwise every candidate substitution matches itself as a modification and
    is dropped from the MTP output.
    """
    mods = gen_mod_dict()
    substitutions = {
        name for rows in mods.values() for name, _, _ in rows if "->" in name
    }
    assert substitutions == set()

    # Concretely: the Ala->Gly delta must have no PTM explanation on alanine.
    delta = aa_subs_ref()["A"]["A to G"]
    assert not [row for row in mods["A"] if row[2] == pytest.approx(delta)]


def test_gen_mod_dict_keeps_real_modifications():
    mods = gen_mod_dict()
    oxidations = [row for row in mods["M"] if row[0] == "Oxidation or Hydroxylation"]
    assert oxidations
    assert oxidations[0][2] == pytest.approx(15.994915, abs=1e-6)


def test_substitution_matrix_is_antisymmetric():
    aas = pd.DataFrame(
        {"one_letter": ["A", "G", "S"], "mono_mass": [71.03711, 57.02146, 87.03203]}
    )
    matrix = calculate_aa_substitution_matrix(aas)
    assert matrix.loc["A", "G"] == pytest.approx(-matrix.loc["G", "A"])
    assert matrix.loc["A", "A"] == 0.0


# -------------------------------------------------------------- column_mapping


def test_column_mapping_is_case_insensitive_and_reorders():
    df = pd.DataFrame({"raw file": ["r"], "CHARGE": [2], "unused": [0]})
    mapped = column_mapping(df, ["Charge", "Raw file"])
    assert list(mapped.columns) == ["Charge", "Raw file"]
    assert mapped.loc[0, "Raw file"] == "r"


def test_column_mapping_skips_absent_columns():
    df = pd.DataFrame({"Charge": [2]})
    assert list(column_mapping(df, ["Charge", "Missing"]).columns) == ["Charge"]


# ----------------------------------------------------------------- ptm/mtp split


def test_ptm_mtp_output_splits_the_candidates(tmp_path):
    (tmp_path / "PTM").mkdir()
    (tmp_path / "MTP").mkdir()

    frame = pd.DataFrame(
        {
            "PTM": ["Oxidation", None, None],
            "mistranslated sequence": ["AAAK", "AAGK", "AAVK"],
            **{
                f"{frame_no}-frame genome substring": [False, False, False]
                for frame_no in range(1, 7)
            },
        }
    )
    # The third candidate is explained by the genome, so it is not an MTP.
    frame.loc[2, "3-frame genome substring"] = True

    ptm_mtp_output(frame, "sample_a", tmp_path)

    ptm = read_frame(tmp_path / "PTM" / "sample_a_PTM")
    mtp = read_frame(tmp_path / "MTP" / "sample_a_MTP")
    assert ptm["mistranslated sequence"].tolist() == ["AAAK"]
    assert mtp["mistranslated sequence"].tolist() == ["AAGK"]


# ------------------------------------------------------------------------ Stage


def test_null_queue_discards_messages():
    assert NullQueue().put(("stdout", "anything")) is None


def test_stage_reads_the_shared_parameters(aas_params):
    stage = Stage(aas_params)
    assert stage.workflow == "MaxQuant"
    assert stage.label_setup == "Label-Free"
    assert stage.label_plex == 1
    assert stage.output_dir.name == "out"


def test_stage_announces_itself_on_the_queue(aas_params):
    queue = RecordingQueue()
    Stage(aas_params, queue)
    assert queue.messages == [("stdout", "Stage initialized.")]


def test_stage_metadata_maps_tmt_channels(aas_params):
    stage = Stage(aas_params)
    assert stage.metadata["MQ"].tolist() == [1, 2]
    assert list(stage.samples) == ["sample_a", "sample_b"]


def test_locate_sample_dir_finds_directories(aas_params):
    stage = Stage(aas_params)
    assert stage._locate_sample_dir("sample_a").name == "sample_a"
    assert stage._locate_sample_dir("sample_a", suffix="_val").name == "sample_a_val"
    assert stage._locate_sample_dir("sample_b") is None


def test_stage_process_sample_is_abstract(aas_params):
    with pytest.raises(NotImplementedError):
        Stage(aas_params).process_sample("sample_a")


def test_record_run_logs_parameters_and_version(aas_params):
    stage = Stage(aas_params)

    log = stage.record_run()

    assert log.name == PROVENANCE_FILE
    entry = json.loads(log.read_text().splitlines()[0])
    assert entry["stage"] == "Stage"
    assert entry["params"]["Detection"]["Protease"] == "Trypsin"
    assert entry["proteolyzer_version"]
    assert entry["timestamp"].endswith("+00:00")


def test_record_run_appends_one_line_per_run(aas_params):
    stage = Stage(aas_params)
    stage.record_run()
    log = stage.record_run()

    assert len(log.read_text().splitlines()) == 2


def test_record_run_is_called_when_a_stage_runs(aas_params):
    from proteolyzer.aas.validation import Validation

    validation = Validation(aas_params)
    validation.run()

    log = Path(aas_params["Utils"]["Output Folder"]) / PROVENANCE_FILE
    assert "Validation" in log.read_text()


def test_detection_dispatches_to_the_workflow_class(aas_params):
    pytest.importorskip("ahocorasick")
    from proteolyzer.aas.detection import Detection, MaxQuant

    detection = Detection(aas_params)
    assert isinstance(detection.detection_workflow, MaxQuant)
    assert detection.detection_workflow.aa_sub_ppm == 10.0
    assert detection.protease == "Trypsin"


def test_detection_rejects_an_unknown_workflow(aas_params):
    pytest.importorskip("ahocorasick")
    from proteolyzer.aas.detection import Detection

    params = dict(aas_params)
    detection = Detection(params)
    detection.workflow = "Spectronaut"
    with pytest.raises(NotImplementedError, match="Spectronaut"):
        detection._initialize_workflow()


def test_validation_reports_missing_inputs(aas_params):
    from proteolyzer.aas.validation import Validation

    queue = RecordingQueue()
    validation = Validation(aas_params, queue)
    validation.process_sample("sample_a")
    assert ("stderr", "Missing data files for sample sample_a") in queue.messages


def test_validation_locates_fragments_spanning_the_substitution(aas_params):
    from proteolyzer.aas.validation import Validation

    validation = Validation(aas_params)
    b_ions, y_ions = validation.frags_containing_aas(
        pd.Series({"mistranslated sequence": "AAGAK", "mistranslated aas positions": 2})
    )
    assert b_ions.tolist() == [3, 4]
    assert y_ions.tolist() == [3, 4]


def test_quantification_reports_missing_inputs(aas_params):
    from proteolyzer.aas.quantification import Quantification

    queue = RecordingQueue()
    Quantification(aas_params, queue).process_sample("sample_a")
    assert ("stderr", "Missing files for sample sample_a") in queue.messages


def test_quantification_label_free_ratio_is_log2(aas_params):
    from proteolyzer.aas.quantification import Quantification

    mtp = pd.DataFrame(
        {
            "DP Base Sequence": ["AAAK"],
            "mistranslated sequence": ["AAGK"],
            "aa subs": ["A:G"],
        }
    )
    mtp_df = pd.DataFrame({"Intensity": [200.0]})
    bp_df = pd.DataFrame({"Intensity": [100.0]})

    quant = Quantification(aas_params)._raas(
        mtp, mtp_df, bp_df, pd.DataFrame({"MQ": [np.nan]}), "Label-Free"
    )
    assert quant["Ratio"].tolist() == [1.0]


def test_quantification_rejects_unknown_label_designation(aas_params):
    from proteolyzer.aas.quantification import Quantification

    with pytest.raises(ValueError, match="Unknown label designation"):
        Quantification(aas_params)._raas(
            pd.DataFrame(
                {
                    "DP Base Sequence": [],
                    "mistranslated sequence": [],
                    "aa subs": [],
                }
            ),
            pd.DataFrame({"Intensity": []}),
            pd.DataFrame({"Intensity": []}),
            pd.DataFrame({"MQ": []}),
            "SILAC",
        )
