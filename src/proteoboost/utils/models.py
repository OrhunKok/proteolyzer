import pandas as pd
from typing import List, Union, Set
from pathlib import Path
import datetime
import logging
import re
from pydantic import BaseModel, Field, field_validator, computed_field
from functools import cached_property
from proteoboost.utils import constants as Constants


class Data(BaseModel):
    file_path: Path = Field(..., description = 'The path to the file to be loaded.')
    load_all_columns: bool = Field(False, description = 'Whether to load all columns from the file.')
    extra_cols_to_load: Union[str, List[str], Set[str]] = Field(None, description = 'Additional columns to load.')

    @field_validator('extra_cols_to_load', mode='before')
    @classmethod
    def _validate_extra_cols_to_load(cls, value: Union[str, List[str], Set[str]]) -> set:
        if isinstance(value, (str)):
            return set([value])
        elif isinstance(value, list) and all(isinstance(item, str) for item in value):
            return set(value) or set()
        elif isinstance(value, set) and all(isinstance(item, str) for item in value):
            return set(value) or set()
        raise TypeError(f"Invalid input type: {type(value)}")
    
    @field_validator('file_path', mode='after')
    @classmethod
    def _validate_file_path(cls, value: Path) -> Path:
        if value.exists():
            return value
        else:
            raise ValueError(f"Path does not exist!")
    
    @computed_field
    @cached_property
    def file_name(self) -> str:
        return self.file_path.stem

    @computed_field
    @cached_property
    def file_extension(self) -> str:
        return self.file_path.suffix
    
    @computed_field
    @cached_property
    def input_type(self) -> str:

        is_diann = self.file_name in Constants.DIANN_FILES and self.file_extension in Constants.DIANN_EXTENSIONS
        is_maxquant = self.file_name in Constants.MAXQUANT_FILES and self.file_extension in Constants.MAXQUANT_EXTENSIONS

        if is_diann and is_maxquant:
            raise ValueError(f"File {self.file_name} with extension {self.file_extension} matches multiple categories.")
        else:
            if is_diann:
                logging.info(f"{self.file_name} determined to be DIA-NN output")
                return 'DIANN'
            elif is_maxquant:
                logging.info(f"{self.file_name} determined to be MaxQuant output")
                return 'MaxQuant'
            else:
                logging.warning(f"{self.file_name} source program could not be determined, certain optimizations will not be performed.")
                return 'Unknown'
            
    @computed_field
    @cached_property
    def cols_subset(self) -> dict:

        if self.input_type != 'Unknown':
            if not self.load_all_columns:
                cols = Constants.SUPPORTED_FILES_COLS_SUBSET[self.input_type][self.file_name]
                if self.extra_cols_to_load:
                    cols = cols | self.extra_cols_to_load
                return cols
            
    @computed_field
    @cached_property
    def cols_rename_mapping(self) -> dict:
        if self.input_type in list(Constants.COLS_RENAME_MAPPING.keys()):
            return Constants.COLS_RENAME_MAPPING[self.input_type]
        
    @computed_field
    @cached_property
    def file_stats(self) -> dict:
        file_stats = self.file_path.stat()
        file_stats = { 
                        'Size (Bytes)': file_stats.st_size,
                        'Created': datetime.datetime.fromtimestamp(file_stats.st_ctime, tz = datetime.timezone.utc).strftime('%Y-%m-%d %H:%M:%S'),
                        'Last Modified': datetime.datetime.fromtimestamp(file_stats.st_mtime, tz = datetime.timezone.utc).strftime('%Y-%m-%d %H:%M:%S'),
                        'Last Accessed': datetime.datetime.fromtimestamp(file_stats.st_atime, tz = datetime.timezone.utc).strftime('%Y-%m-%d %H:%M:%S')
                     }
        return file_stats


class ProcessedData(BaseModel):
    data: pd.DataFrame = Field(..., description='The processed DataFrame.')
    ID_COL : str = Field(..., description='Reference used for unique IDs')
    LABEL_FREE : bool = Field(..., description='Data is label-free.')
    LABEL_GROUP_CAPTURE : re.Pattern = Field(..., description='Regex pattern used to determine labelling.')

    class Config:
        arbitrary_types_allowed = True  # Allow Pandas DataFrames

    def __init__(self, processor: "DataProcessor", **kwargs):

        processor_attrs = {k: getattr(processor, k) for k in processor.__slots__ if k in self.__fields__}
        super().__init__(**processor_attrs)

    @computed_field
    @cached_property
    def unique_runs(self) -> set:
        if 'Run' in self.data.columns:
            return set(self.data['Run'].unique())
        else:
            return set()

    @computed_field
    @cached_property
    def unique_ids(self) -> int:
        return self.data[self.ID_COL].nunique()
