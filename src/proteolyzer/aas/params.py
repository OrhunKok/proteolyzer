"""Schema and loader for the AAS pipeline parameter file.

The YAML file is validated against :class:`ParamsSchema` and then flattened:
each section keeps its shared keys and is merged with the block named after
``Utils.Workflow``, so a stage can read ``params["Detection"]["ALT ppm"]``
without knowing which search engine produced the data.
"""

from pathlib import Path
from typing import Any, Literal

import yaml
from pydantic import AliasChoices, BaseModel, Field


class ParamsModel(BaseModel):
    model_config = {
        "populate_by_name": True,
        "extra": "forbid",
    }


class Utils(ParamsModel):
    Data_Folder: Path = Field(alias="Data Folder")
    Output_Folder: Path = Field(alias="Output Folder")
    Metadata_File: Path = Field(alias="Metadata File")
    Workflow: Literal["MaxQuant"]
    Labelling_Setup: Literal["Label-Free", "TMT"] = Field(alias="Labelling Setup")
    Label_Plex: int = Field(alias="Label Plex")


class Translation(ParamsModel):
    Genome_FASTA: Path = Field(alias="Genome FASTA")
    Translated_Frames_Folder: Path = Field(alias="Translated Frames Folder")


class DetectionMaxQuant(ParamsModel):
    Detection_PEP: float = Field(alias="Detection PEP")
    AA_Substitution_ppm: float = Field(alias="AA Substitution ppm")
    # "PTM ppm" is the old name; ALT is the peptide class it assigns.
    ALT_ppm: float = Field(
        validation_alias=AliasChoices("ALT ppm", "PTM ppm"),
        serialization_alias="ALT ppm",
    )
    Positional_Probability_Threshold: float = Field(
        alias="Positional Probability Threshold"
    )
    C_n_term_Modification_Threshold: float = Field(
        alias="C/n-term Modification Threshold"
    )


class Detection(ParamsModel):
    Protease: str
    Protein_FASTA: Path = Field(alias="Protein FASTA")
    MaxQuant: DetectionMaxQuant


class ValidationMaxQuant(ParamsModel):
    Validation_PEP: float = Field(alias="Validation PEP")
    PIF: float
    Fragment_Evidence: float = Field(alias="Fragment Evidence")


class Validation(ParamsModel):
    MaxQuant: ValidationMaxQuant


class QuantificationMaxQuant(ParamsModel):
    Minimum_Quantity: int = Field(alias="Minimum Quantity")


class Quantification(ParamsModel):
    MaxQuant: QuantificationMaxQuant


class ParamsSchema(ParamsModel):
    Utils: Utils
    Translation: Translation
    Detection: Detection
    Validation: Validation
    Quantification: Quantification


def _load_yaml(filepath: str | Path) -> dict[str, Any]:
    """Load raw parameters from a YAML file using safe_load."""
    with open(filepath) as f:
        return yaml.safe_load(f)


def load_params(params: str | Path | dict) -> dict:
    """Validate parameters and merge the workflow-specific sections.

    Accepts a path to a YAML file or an already-loaded mapping.
    """
    if isinstance(params, dict):
        loaded_params = params
    elif isinstance(params, (str, Path)):
        if not Path(params).is_file():
            raise ValueError(f"The file {params} does not exist.")
        loaded_params = _load_yaml(params)
    else:
        raise TypeError(
            f"params must be a file path or dict, got {type(params).__name__}"
        )

    validated_params = ParamsSchema.model_validate(loaded_params).model_dump(
        by_alias=True
    )
    workflow = validated_params["Utils"]["Workflow"]

    merged_params = {}
    for section_name, section_data in validated_params.items():
        if isinstance(section_data, dict):
            base = {k: v for k, v in section_data.items() if not isinstance(v, dict)}
            workflow_section = section_data.get(workflow, {})
            merged_params[section_name] = {**base, **workflow_section}
        else:
            merged_params[section_name] = section_data

    return merged_params
