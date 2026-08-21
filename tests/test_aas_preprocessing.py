"""The preprocessor turns MaxQuant text exports into the pipeline's parquet inputs."""

from pathlib import Path

import pandas as pd
import pytest

pytest.importorskip("yaml")
pytest.importorskip("fastparquet")

from proteolyzer.aas.preprocessing import Preprocessor  # noqa: E402


@pytest.fixture
def search_dir(aas_params):
    """A MaxQuant "Detection" search directory (it has dependentPeptides.txt)."""
    search = Path(aas_params["Utils"]["Data Folder"]) / "sample_a"
    (search / "dependentPeptides.txt").write_text("marker\n")

    pd.DataFrame(
        {
            "Raw file": ["r1", "r2", "r3"],
            "Charge": [2, 2, 2],
            "m/z": [500.0, 600.0, 700.0],
            "Mass": [1000.0, 1200.0, 1400.0],
            "Retention time": [10.0, 11.0, 12.0],
            "Intensity": [100, 200, 300],
            "DP Mass Difference": [14.0, None, 28.0],
            "DP PEP": [0.001, 0.001, 0.5],
            "DP Decoy": [None, None, None],
            "DP Proteins": ["P1;P2", "P3", "CON__P4"],
            "Reverse": [None, None, None],
            "DP Base Sequence": ["AAAK", "BBBK", "CCCK"],
            "DP Probabilities": ["AA(0.9)AK", "BB(0.9)BK", "CC(0.9)CK"],
            "DP Positional Probability": [0.9, 0.9, 0.9],
            "DP Base Scan Number": [1, 2, 3],
            "DP Mod Scan Number": [1, 2, 3],
            "MSMS Scan Numbers": ["1;2", "3", "4"],
        }
    ).to_csv(search / "allPeptides.txt", sep="\t", index=False)

    pd.DataFrame(
        {
            "Raw file": ["r1", "r2"],
            "Charge": [2, 2],
            "m/z": [500.0, 600.0],
            "Mass": [1000.0, 1200.0],
            "Retention time": [10.0, 11.0],
            "Reverse": [None, "+"],
            "Potential contaminant": [None, None],
            "Sequence": ["AAAK", "BBBK"],
            "PIF": [0.9, 0.9],
            "PEP": [0.001, 0.001],
            "Mass error [ppm]": [1.0, 2.0],
            "MS/MS scan number": [1, 3],
            "Intensity": [100, 200],
        }
    ).to_csv(search / "evidence.txt", sep="\t", index=False)

    pd.DataFrame(
        {
            "Sequence": ["AAAK", "BBBK", "CCCK"],
            "Start position": [1, 5, 9],
            "Amino acid after": ["G", "-", "G"],
            "Amino acid before": ["-", "K", "K"],
            "Reverse": [None, None, "+"],
            "Potential contaminant": [None, None, None],
        }
    ).to_csv(search / "peptides.txt", sep="\t", index=False)

    return search


def test_run_converts_the_needed_files_only(aas_params, search_dir):
    Preprocessor.MaxQuant(aas_params).run()

    converted = {p.name for p in search_dir.glob("*.parquet")}
    assert converted == {"allPeptides.parquet", "evidence.parquet", "peptides.parquet"}


def test_run_creates_the_output_layout(aas_params, search_dir):
    Preprocessor.MaxQuant(aas_params).run()

    output_dir = Preprocessor.MaxQuant(aas_params).output_dir
    assert (output_dir / "ALT").is_dir()
    assert (output_dir / "SAAP").is_dir()


def test_allpeptides_filtering(aas_params, search_dir):
    Preprocessor.MaxQuant(aas_params).run()
    frame = pd.read_parquet(search_dir / "allPeptides.parquet", engine="fastparquet")

    # r2 has no mass difference, r3 is a contaminant above the PEP threshold.
    assert frame["Raw file"].tolist() == ["r1"]
    assert frame.loc[0, "MS/MS scan number"] == [1, 2]
    assert frame.loc[0, "Leading.Razor.DP.Protein"] == "P1"


def test_evidence_drops_reverse_hits(aas_params, search_dir):
    Preprocessor.MaxQuant(aas_params).run()
    frame = pd.read_parquet(search_dir / "evidence.parquet", engine="fastparquet")
    assert frame["Sequence"].tolist() == ["AAAK"]


def test_peptides_are_labelled_by_terminus(aas_params, search_dir):
    Preprocessor.MaxQuant(aas_params).run()
    frame = pd.read_parquet(search_dir / "peptides.parquet", engine="fastparquet")
    assert frame.set_index("Sequence")["Terminus"].to_dict() == {
        "AAAK": "N",
        "BBBK": "C",
    }


def test_missing_data_folder_is_reported(aas_params, tmp_path):
    """Regression: this surfaced as a bare errno 2 from iterdir()."""
    aas_params["Utils"]["Data Folder"] = str(tmp_path / "not_there")

    with pytest.raises(FileNotFoundError, match="Data Folder does not exist"):
        Preprocessor.MaxQuant(aas_params).run()


def test_empty_data_folder_is_reported(aas_params, tmp_path):
    empty = tmp_path / "empty"
    empty.mkdir()
    aas_params["Utils"]["Data Folder"] = str(empty)

    with pytest.raises(FileNotFoundError, match="No MaxQuant search output"):
        Preprocessor.MaxQuant(aas_params).run()
