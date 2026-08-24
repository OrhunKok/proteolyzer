"""Engine format quirks moved from streamlit-DO-MS.

These read a column exactly as the engine writes it, ahead of proteolyzer's
own rename -- so the columns below use the engine's own names, not the
canonical schema.
"""

import pandas as pd

from proteolyzer.core.quirks import (
    contaminant_flag,
    decoy_flag,
    fragpipe_run_names,
    psm_identified,
    run_name_from_path,
)

# --------------------------------------------------------------- run_name_from_path


def test_run_name_from_path_takes_the_basename():
    assert run_name_from_path("/data/runs/dOK083.raw") == "dOK083.raw"


def test_run_name_from_path_handles_windows_separators():
    assert run_name_from_path(r"C:\data\runs\dOK083.raw") == "dOK083.raw"


def test_run_name_from_path_strips_a_suffix():
    assert (
        run_name_from_path("interact-dOK083.pepXML", suffix=".pepXML")
        == "interact-dOK083"
    )


def test_run_name_from_path_strips_a_prefix():
    assert run_name_from_path("interact-dOK083", prefix="interact-") == "dOK083"


def test_run_name_from_path_strips_prefix_and_suffix_together():
    path = "/data/interact-dOK083.pepXML"
    assert run_name_from_path(path, prefix="interact-", suffix=".pepXML") == "dOK083"


def test_run_name_from_path_leaves_an_undecorated_name_alone():
    assert (
        run_name_from_path("dOK083", prefix="interact-", suffix=".pepXML") == "dOK083"
    )


# --------------------------------------------------------------- fragpipe_run_names


def test_fragpipe_run_names_takes_everything_before_the_last_three_fields():
    spectrum = pd.Series(["dOK083.d.1234.5678.2"])
    spectrum_file = pd.Series(["interact-dOK083.pepXML"])

    assert fragpipe_run_names(spectrum, spectrum_file).tolist() == ["dOK083.d"]


def test_fragpipe_run_names_falls_back_to_spectrum_file():
    spectrum = pd.Series([None])
    spectrum_file = pd.Series(["interact-dOK084.pepXML"])

    assert fragpipe_run_names(spectrum, spectrum_file).tolist() == ["dOK084"]


def test_fragpipe_run_names_mixes_rows_with_and_without_spectrum():
    spectrum = pd.Series(["dOK083.1.1.2", None])
    spectrum_file = pd.Series(["interact-dOK083.pepXML", "interact-dOK084.pepXML"])

    assert fragpipe_run_names(spectrum, spectrum_file).tolist() == [
        "dOK083",
        "dOK084",
    ]


# --------------------------------------------------------------- psm_identified


def test_psm_identified_is_true_for_a_real_sequence():
    assert psm_identified(pd.Series(["PEPTIDER"])).tolist() == [True]


def test_psm_identified_is_false_for_a_single_space():
    assert psm_identified(pd.Series([" "])).tolist() == [False]


def test_psm_identified_is_false_for_an_empty_string():
    assert psm_identified(pd.Series([""])).tolist() == [False]


def test_psm_identified_is_false_for_null():
    assert psm_identified(pd.Series([None])).tolist() == [False]


def test_psm_identified_reads_all_three_the_same_way():
    result = psm_identified(pd.Series(["PEPTIDER", " ", "", None]))
    assert result.tolist() == [True, False, False, False]


# --------------------------------------------------------------- decoy_flag


def test_decoy_flag_of_real_bools():
    assert decoy_flag(pd.Series([True, False])).tolist() == [True, False]


def test_decoy_flag_of_numeric_zero_and_one():
    assert decoy_flag(pd.Series([1, 0])).tolist() == [True, False]


def test_decoy_flag_of_the_word():
    result = decoy_flag(pd.Series(["True", "False", "TRUE", "false"]))
    assert result.tolist() == [True, False, True, False]


def test_decoy_flag_of_numeric_text():
    assert decoy_flag(pd.Series(["1", "0"])).tolist() == [True, False]


def test_decoy_flag_of_a_categorical_column():
    column = pd.Series(["True", "False", "True"], dtype="category")
    assert decoy_flag(column).tolist() == [True, False, True]


# --------------------------------------------------------------- contaminant_flag


def test_contaminant_flag_of_a_contaminant():
    protein = pd.Series(["contam_sp|P02768|ALBU_HUMAN"])
    assert contaminant_flag(protein).tolist() == [True]


def test_contaminant_flag_of_a_real_protein():
    protein = pd.Series(["sp|P02768|ALBU_HUMAN"])
    assert contaminant_flag(protein).tolist() == [False]


def test_contaminant_flag_accepts_a_different_mark():
    protein = pd.Series(["CON__P02768", "sp|P02768"])
    assert contaminant_flag(protein, mark="CON__").tolist() == [True, False]
