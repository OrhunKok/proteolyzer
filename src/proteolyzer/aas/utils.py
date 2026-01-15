import os
import pandas as pd
import numpy as np
import pickle
import yaml
from pydantic import BaseModel, Field
from typing import Dict, Union, Any, Literal
from multiprocessing import Queue
from importlib import resources
from pathlib import Path
from .config import Config

CONFIG = Config()

class NullQueue:
    """A 'do-nothing' queue to replace the real Queue when multiprocessing isn't used."""
    def put(self, item):
        pass

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
    Labelling_Setup: Literal["Label-Free", 'TMT'] = Field(alias="Labelling Setup")
    Label_Plex: int = Field(alias="Label Plex")

class Translation(ParamsModel):
    Genome_FASTA: Path = Field(alias="Genome FASTA")
    Translated_Frames_Folder: Path = Field(alias="Translated Frames Folder")

class DetectionMaxQuant(ParamsModel):
    Detection_PEP: float = Field(alias="Detection PEP")
    AA_Substitution_ppm: float = Field(alias="AA Substitution ppm")
    PTM_ppm: float = Field(alias="PTM ppm")
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

def _load_yaml(filepath: str) -> Dict[str, Any]:
    """Load raw parameters from a YAML file using safe_load."""
    with open(filepath, 'r') as f:
        return yaml.safe_load(f)

def _load_params(params: Union[str, Dict]) -> dict:
    """Load YAML or dict params and merge workflow-specific sections automatically."""

    if isinstance(params, dict):
        loaded_params = params
    elif os.path.isfile(params):
        loaded_params = _load_yaml(params)
    else:
        raise ValueError(f"The file {params} does not exist.")
    
    validated_params = ParamsSchema.model_validate(loaded_params).model_dump(by_alias=True)
    workflow = validated_params['Utils']['Workflow']
    
    merged_params = {}
    for section_name, section_data in validated_params.items():
        if isinstance(section_data, dict):
            base = {k: v for k, v in section_data.items() if not isinstance(v, dict)}
            workflow_section = section_data.get(workflow, {})
            merged_params[section_name] = {**base, **workflow_section}
        else:
            merged_params[section_name] = section_data

    return merged_params


def column_mapping(df: pd.DataFrame, cols2keep: list) -> pd.DataFrame:
    """Standardized column selection + renaming based on 'cols2keep'."""
    col_lookup = {col.lower(): col for col in df.columns}
    rename_map = {
        col_lookup[col.lower()]: col
        for col in cols2keep
        if col.lower() in col_lookup
    }

    df = df[list(rename_map.keys())].rename(columns=rename_map)
    df = df[[col for col in cols2keep if col in df.columns]]
    return df


def aa_subs_ref(
               aa_subs_file : str = r"unimod_modifications.csv"
               ) -> dict:
    with resources.files("proteolyzer.resources").joinpath(aa_subs_file).open("r") as f:
            aa_subs = pd.read_csv(f)

    aa_subs = aa_subs[aa_subs['classification'] == 'AA substitution']
    aa_subs['code_name'] = aa_subs['code_name'].str.replace('2', '->')
    aa_subs['sub_aa'] = aa_subs['code_name'].str.split(' ').str[0].str.split('->').str[1]
    aa_subs['sub_aa'] = aa_subs['sub_aa'].map(CONFIG.AminoAcids.CODE_TO_SYMBOL)
    aa_subs['compact_name'] = aa_subs['one_letter'] + ' to ' + aa_subs['sub_aa']
    subs_ref = aa_subs[['one_letter', 'compact_name', 'mono_mass']].groupby('one_letter').apply(lambda x: x.set_index('compact_name')['mono_mass'].to_dict()).to_dict()
    return subs_ref


def gen_mod_dict(
                mod_file : str = r"unimod_modifications.csv"
                ) -> dict:
    with resources.files("proteolyzer.resources").joinpath(mod_file).open("r") as f:
        mod_df = pd.read_csv(f, sep = ',', usecols = ['full_name', 'position', 'one_letter', 'mono_mass'])

    all_mods_dict = {}
    for k, g in mod_df.groupby('one_letter'):
        all_mods_dict[k] = g.drop('one_letter', axis = 1)[['full_name', 'position', 'mono_mass']].values.tolist()

    return all_mods_dict

class Preprocessor:

    class MaxQuant:
        FILES            = CONFIG.MaxQuant.FILES
        FILES_NEEDED     = CONFIG.MaxQuant.FILES_NEEDED
        LOAD_COLS        = CONFIG.MaxQuant.LOAD_COLS
        RENAME_MAP       = CONFIG.MaxQuant.COLS_RENAME_MAPPING

        def __init__(self, params: Union[str, Dict], queue: Queue = None):
            if isinstance(params, str):
                self.params = _load_params(params)
            elif isinstance(params, dict):
                self.params = params
            else:
                raise ValueError("params must be a file path or dict")
            
            self.queue = queue if queue is not None else NullQueue()

            self.data_folder = self.params.get('Utils').get('Data Folder')
            self.translated_frames = self.params.get('Translation').get('Translated Frames Folder')
            self.output_dir = self.params.get('Utils').get('Output Folder')
            self.label_setup = self.params.get('Utils').get('Labelling Setup')
            self.label_plex = self.params.get('Utils').get('Label Plex')

        def run(self):

            self.translated_frames.mkdir(parents=True, exist_ok=True)
            (self.output_dir / 'PTM').mkdir(parents=True, exist_ok=True)
            (self.output_dir / 'MTP').mkdir(parents=True, exist_ok=True)

            for root, dirs, _ in os.walk(self.data_folder):
                for dirname in dirs:
                    subdir = Path(root) / dirname
                    files = {f for f in os.listdir(subdir)}

                    # Detect type of search
                    search_type = "Detection" if "dependentPeptides.txt" in files else "Validation"
                    needed = {f + ".txt" for f in self.FILES_NEEDED[search_type]}

                    for filename in files:
                        if filename not in needed:
                            continue

                        txt_path = subdir / filename
                        base = filename.replace(".txt", "")
                        parquet_path = subdir / f"{base}.parquet"

                        cols2keep = self.LOAD_COLS[base]
                        df = pd.read_csv(txt_path, delimiter="\t")

                        if (
                            base == "evidence"
                            and self.label_setup == "TMT"
                            and "_val" in dirname
                        ):
                            cols2keep = cols2keep + [
                                f"Reporter intensity corrected {i}"
                                for i in range(1, self.label_plex + 1)
                            ]

                        df = column_mapping(df, cols2keep)

                        if base == "allPeptides":
                            df = self._allpeptides(df)

                        elif base == "evidence":
                            df = self._evidence(df)

                        elif base == "msms":
                            df = self._msms(df)

                        elif base == "peptides":
                            df = self._peptides(df)

                        df.reset_index(drop=True, inplace=True)
                        df.to_parquet(parquet_path, engine="fastparquet")
                        print(f"{dirname} Converted {txt_path} → {parquet_path}")


        def _evidence(self, df: pd.DataFrame):
            df = df[(df["Reverse"].isna()) & (df["Potential contaminant"].isna())]
            return df

        def _allpeptides(self, df: pd.DataFrame):
            df["Leading.Razor.DP.Protein"] = df["DP Proteins"].str.split(";").str[0]
            df["Con.DP.Protein"] = df["DP Proteins"].str.contains("CON__", na=False)

            df = df[(df["Reverse"].isna()) & (~df["Con.DP.Protein"])]

            df = df.loc[
                ~np.isnan(df["DP Mass Difference"])
                & (df["DP PEP"] <= CONFIG.PEP_THRESHOLD)
            ]

            df = df[~df["MSMS Scan Numbers"].isna()]
            df["MSMS Scan Numbers"] = df["MSMS Scan Numbers"].str.split(";").apply(
                lambda x: list(map(int, x))
            )

            df.rename(columns={"MSMS Scan Numbers": "MS/MS scan number"}, inplace=True)
            return df

        def _msms(self, df: pd.DataFrame):
            df = df[df["Reverse"].isna()]
            df.rename(columns={"Scan number": "MS/MS scan number"}, inplace=True)

            df["Matches"] = df["Matches"].str.split(";")
            df = df[df["Matches"].notna()]

            df["Matches"] = df.apply(
                lambda row: [
                    frag for frag in row["Matches"]
                    if "a" not in frag and "(" not in frag
                    and "NH3" not in frag and "H2O" not in frag
                ],
                axis=1,
            )

            df["Frag.Type"] = df["Matches"].apply(lambda m: [frag[0] for frag in m])
            df["Frag.Number"] = df["Matches"].apply(lambda m: [int(frag[1:]) for frag in m])

            df = df.drop(columns="Matches").explode(["Frag.Type", "Frag.Number"])
            return df

        def _peptides(self, df: pd.DataFrame):
            df = df[(df["Reverse"].isna()) & (df["Potential contaminant"].isna())]
            df = df[(df["Start position"] == 1) | (df["Amino acid after"] == "-")]
            
            df['Terminus'] = np.where(df["Amino acid after"] == "-", "C", "N")
            return df


    class DIANN:
        FILES = CONFIG.DIANN.FILES
        FILE_EXT = CONFIG.DIANN.FILE_EXTENSIONS

        def __init__(self, data_folder: str):
            self.data_folder = Path(data_folder)

        def run(self):
            print("[DIANN] No preprocessing pipeline implemented yet.")
            pass


def ptm_mtp_output(
                  dp_df : pd.DataFrame,
                  sample_name : str,
                  output_dir : Path,
                  ):
    
    ptm_df = dp_df.copy()
    ptm_df = ptm_df[ptm_df['PTM'].notna()] 
    with open(output_dir / 'PTM' / f'{sample_name}_PTM.p', 'wb') as f:
        pickle.dump(ptm_df, f)

    # MTP are peptides with potential AAS that cannot be explained by PTM and do not have homologous sequence in genome
    mtp_df = dp_df.copy()
    mtp_df = mtp_df[mtp_df['PTM'].isna()] 
    mtp_df = mtp_df[~mtp_df['mistranslated sequence'].isnull()]
    mtp_df = mtp_df[~mtp_df.apply(lambda x: any([x['1-frame genome substring'], x['2-frame genome substring'], x['3-frame genome substring'],
                                                x['4-frame genome substring'], x['5-frame genome substring'], x['6-frame genome substring']]), axis = 1)]
    with open(output_dir / 'MTP' / f'{sample_name}_MTP.p', 'wb') as f:
        pickle.dump(mtp_df, f)


def calculate_aa_substitution_matrix(processed_amino_acids_df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculates the pairwise mass difference matrix (Row AA mass - Column AA mass).
    """
    aa_vector: pd.Series = processed_amino_acids_df.set_index('one_letter')['mono_mass']
    
    aa_subs_pairwise: pd.DataFrame = pd.DataFrame(
        aa_vector.values[np.newaxis, :] - aa_vector.values[:, np.newaxis],
        index=aa_vector.index.values,
        columns=aa_vector.index.values
    )

    return aa_subs_pairwise