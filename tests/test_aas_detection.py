"""End-to-end smoke test for the detection stage.

Builds one synthetic dependent peptide whose mass shift matches an Ala->Gly
substitution, runs the stage that classifies it, and checks the artefacts it
writes. It is deliberately structural: it guards the plumbing across
detection.py, not the numeric thresholds of the method.
"""

import pickle
from pathlib import Path

import pandas as pd
import pytest

pytest.importorskip("yaml")
pytest.importorskip("fastparquet")
pytest.importorskip("ahocorasick")

from proteolyzer.aas.detection import Detection  # noqa: E402
from proteolyzer.aas.io import read_frame  # noqa: E402
from proteolyzer.aas.preprocessing import Preprocessor  # noqa: E402
from proteolyzer.aas.utils import aa_subs_ref  # noqa: E402

#: Mass shift of the substitution the fixture encodes, and the peptides it maps.
ALA_TO_GLY = aa_subs_ref()["A"]["A to G"]
BASE_SEQUENCE = "AAAK"
MISTRANSLATED = "AGAK"


@pytest.fixture
def detection_inputs(aas_params):
    """A search directory plus translated frames, as the earlier stages leave them."""
    search = Path(aas_params["Utils"]["Data Folder"]) / "sample_a"
    (search / "dependentPeptides.txt").write_text("marker\n")

    pd.DataFrame(
        {
            "Raw file": ["r1"],
            "Charge": [2],
            "m/z": [500.0],
            "Mass": [1000.0],
            "Retention time": [10.0],
            "Intensity": [1000],
            # Matches the Ala->Gly delta exactly, so it is inside any tolerance.
            "DP Mass Difference": [ALA_TO_GLY],
            "DP PEP": [0.001],
            "DP Decoy": [None],
            "DP Proteins": ["P00001"],
            "Reverse": [None],
            "DP Base Sequence": [BASE_SEQUENCE],
            # The localization probability sits above the configured threshold.
            "DP Probabilities": ["AA(0.95)AK"],
            "DP Positional Probability": [0.95],
            "DP Base Scan Number": [1],
            "DP Mod Scan Number": [1],
            "MSMS Scan Numbers": ["1;2"],
        }
    ).to_csv(search / "allPeptides.txt", sep="\t", index=False)

    # Two rows with different mass errors: the posterior needs a non-zero spread.
    pd.DataFrame(
        {
            "Raw file": ["r1", "r1"],
            "Charge": [2, 2],
            "m/z": [500.0, 500.0],
            "Mass": [1000.0, 1000.0],
            "Retention time": [10.0, 11.0],
            "Reverse": [None, None],
            "Potential contaminant": [None, None],
            "Sequence": [BASE_SEQUENCE, MISTRANSLATED],
            "PIF": [0.9, 0.9],
            "PEP": [0.001, 0.002],
            "Mass error [ppm]": [1.0, 3.0],
            "MS/MS scan number": [1, 2],
            "Intensity": [1000, 100],
        }
    ).to_csv(search / "evidence.txt", sep="\t", index=False)

    pd.DataFrame(
        {
            "Sequence": [BASE_SEQUENCE],
            "Start position": [1],
            "Amino acid after": ["G"],
            "Amino acid before": ["-"],
            "Reverse": [None],
            "Potential contaminant": [None],
        }
    ).to_csv(search / "peptides.txt", sep="\t", index=False)

    Preprocessor.MaxQuant(aas_params).run()

    # Stand in for FrameTranslator: six frames with no homologue of the
    # mistranslated peptide, so the candidate is not explained by the genome.
    frames = Path(aas_params["Translation"]["Translated Frames Folder"])
    for frame in range(1, 7):
        for name in (f"frame_{frame}.p", f"frame_{frame}_il_ambigous.p"):
            with open(frames / name, "wb") as f:
                pickle.dump("MTTTTVVVVWWWW*", f)

    return search


def test_detection_writes_ptm_and_mtp_artefacts(aas_params, detection_inputs):
    detector = Detection(aas_params)
    detector.run()

    output_dir = detector.output_dir
    assert (output_dir / "PTM" / "sample_a_PTM.parquet").exists()
    assert (output_dir / "MTP" / "sample_a_MTP.parquet").exists()
    assert (output_dir / "MTP" / "sample_a_MTP_Filtered_Stage_1.parquet").exists()
    assert (output_dir / "MTP" / "sample_a_FASTA.parquet").exists()
    assert (output_dir / "sample_a_validation.fasta").exists()


def test_detection_identifies_the_substitution(aas_params, detection_inputs):
    detector = Detection(aas_params)
    detector.run()

    mtp = read_frame(detector.output_dir / "MTP" / "sample_a_MTP")

    assert mtp["aa subs"].tolist() == ["A to G"]
    assert mtp["mistranslated sequence"].tolist() == [MISTRANSLATED]
    assert mtp["mistranslated aas positions"].tolist() == [1]
    # No PTM explains the shift, and no frame contains the peptide.
    assert mtp["PTM"].isna().all()
    assert not mtp[[f"{i}-frame genome substring" for i in range(1, 7)]].any().any()


def test_detection_flags_candidates_found_in_the_genome(aas_params, detection_inputs):
    """A homologous sequence in any frame disqualifies the candidate.

    That leaves no candidates at all, which the validation step has to survive.
    """
    frames = Path(aas_params["Translation"]["Translated Frames Folder"])
    with open(frames / "frame_3.p", "wb") as f:
        # Preceded by a cleavage site, as the search requires.
        pickle.dump(f"MTTTK{MISTRANSLATED}VVVV", f)

    detector = Detection(aas_params)
    detector.run()

    assert read_frame(detector.output_dir / "MTP" / "sample_a_MTP").empty
    stage_1 = detector.output_dir / "MTP" / "sample_a_MTP_Filtered_Stage_1"
    assert read_frame(stage_1).empty
    # The validation FASTA is still written, holding just the input proteins.
    assert (
        detector.output_dir / "sample_a_validation.fasta"
    ).read_text() == ">sp|P00001|TEST\nMAIV\n"


def test_write_fasta_appends_candidates_to_the_protein_fasta(
    aas_params, detection_inputs
):
    workflow = Detection(aas_params).detection_workflow

    workflow.write_fasta(
        pd.DataFrame(
            {
                "DP Base Sequence": [BASE_SEQUENCE],
                "mistranslated sequence": [MISTRANSLATED],
                "destination aa": ["G"],
                "mistranslated aas positions": [1],
                "aa subs": ["A to G"],
                "Leading.Razor.DP.Protein": ["P00001"],
            }
        ),
        "fasta_sample",
    )

    fasta = (workflow.output_dir / "fasta_sample_validation.fasta").read_text()
    # The input protein FASTA is copied first, then candidates are appended.
    assert fasta == (
        ">sp|P00001|TEST\nMAIV\n"
        f">MTP|({BASE_SEQUENCE})(1)(A:G)(P00001)\n{MISTRANSLATED}\n"
    )


def test_write_fasta_expands_ambiguous_xle_substitutions(aas_params, detection_inputs):
    """A destination of J (Xle) is emitted as both isoleucine and leucine."""
    workflow = Detection(aas_params).detection_workflow

    workflow.write_fasta(
        pd.DataFrame(
            {
                "DP Base Sequence": ["AAAK"],
                "mistranslated sequence": ["AJAK"],
                "destination aa": ["J"],
                "mistranslated aas positions": [1],
                "aa subs": ["A to J"],
                "Leading.Razor.DP.Protein": ["P00001"],
            }
        ),
        "xle_sample",
    )

    fasta_df = read_frame(workflow.output_dir / "MTP" / "xle_sample_FASTA")
    assert set(fasta_df["destination aa"]) == {"I", "L"}
    assert set(fasta_df["mistranslated sequence"]) == {"AIAK", "ALAK"}
    assert set(fasta_df["aa subs"]) == {"A:I", "A:L"}
