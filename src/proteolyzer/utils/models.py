import pandas as pd
from typing import List, Union, Set
from pathlib import Path
import datetime
import logging
from pydantic import BaseModel, Field, field_validator, computed_field
from typing import Optional
from functools import cached_property, cache
from .config import Config

CONFIG = Config()

class Data(BaseModel):
    file_path: Path = Field(..., description="The path to the file to be loaded.")
    load_all_columns: bool = Field(
        False, description="Whether to load all columns from the file."
    )
    extra_cols_to_load: Union[str, List[str], Set[str]] = Field(
        None, description="Additional columns to load."
    )
    INPUT_TYPE: Optional[str] = Field(
        None, description="Manually set the input data type (e.g., 'DIANN', 'MaxQuant')."
    )

    @field_validator("extra_cols_to_load", mode="before")
    @classmethod
    def _validate_extra_cols_to_load(
        cls, value: Union[str, List[str], Set[str]]
    ) -> set:
        if isinstance(value, str):
            return set([value])
        elif isinstance(value, list) and all(isinstance(item, str) for item in value):
            return set(value) or set()
        elif isinstance(value, set) and all(isinstance(item, str) for item in value):
            return set(value) or set()
        raise TypeError(f"Invalid input type: {type(value)}")

    @field_validator("file_path", mode="after")
    @classmethod
    def _validate_path(cls, value: Path) -> Path:
        if value.exists():
            return value
        else:
            raise ValueError("Path does not exist!")

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
        user_override = self.INPUT_TYPE
        
        is_diann = (self.file_name in CONFIG.DIANN.FILES and self.file_extension in CONFIG.DIANN.FILE_EXTENSIONS)
        is_maxquant = (self.file_name in CONFIG.MaxQuant.FILES and self.file_extension in CONFIG.MaxQuant.FILE_EXTENSIONS)

        if is_diann and is_maxquant:
            raise ValueError(
                f"File {self.file_name} with extension {self.file_extension} matches multiple categories."
            )
        
        if is_diann:
            auto_type = "DIANN"
        elif is_maxquant:
            auto_type = "MaxQuant"
        else:
            auto_type = "Unknown"
        
        if user_override not in (None, "Unknown"):
            if auto_type != "Unknown" and user_override != auto_type:
                logging.warning(
                    f"User input '{user_override}' conflicts with file type '{auto_type}'. Recommend using auto-detected type."
                )
            logging.info(f"Using manually set input type: {user_override}")
            return user_override
            
        if auto_type == "Unknown":
            logging.warning(
                f"{self.file_name} source program could not be determined, certain optimizations will not be performed."
            )
        else:
            logging.info(f"{self.file_name} determined to be {auto_type} output")
            
        return auto_type

    @computed_field
    @cached_property
    def cols_subset(self) -> set:
        if self.input_type == "Unknown" or self.load_all_columns:
            return None
            
        cols = getattr(getattr(CONFIG, self.input_type, None), 'LOAD_COLS', {}).get(self.file_name)

        if cols is None:
            return None
        
        if self.extra_cols_to_load:
            return set(cols) | set(self.extra_cols_to_load)
        
        return cols

    @computed_field
    @cached_property
    def cols_rename_mapping(self) -> dict:
        try:
            config_block = getattr(CONFIG, self.input_type)
            return config_block.COLS_RENAME_MAPPING
        except AttributeError:
            return {}

    @computed_field
    @cached_property
    def file_stats(self) -> dict:
        file_stats = self.file_path.stat()
        file_stats = {
            "Size (Bytes)": file_stats.st_size,
            "Created": datetime.datetime.fromtimestamp(
                file_stats.st_ctime, tz=datetime.timezone.utc
            ).strftime("%Y-%m-%d %H:%M:%S"),
            "Last Modified": datetime.datetime.fromtimestamp(
                file_stats.st_mtime, tz=datetime.timezone.utc
            ).strftime("%Y-%m-%d %H:%M:%S"),
            "Last Accessed": datetime.datetime.fromtimestamp(
                file_stats.st_atime, tz=datetime.timezone.utc
            ).strftime("%Y-%m-%d %H:%M:%S"),
        }
        return file_stats
    
    def load(self):
        from .loader import DataLoader
        return LoadedData(DataLoader(self))


class ProcessedData(pd.DataFrame): 
    """
    A specialized DataFrame to hold processed data and its metadata.
    Inherits all pandas DataFrame methods and attributes.
    """
    _metadata = ['ID_COL', 'LABEL_FREE', 'LABEL_GROUP_CAPTURE', 'PROTEASE']

    @property
    def _constructor(self):
        return ProcessedData

    def __init__(self, data=None, ID_COL=None, LABEL_FREE=None, 
                 LABEL_GROUP_CAPTURE=None, PROTEASE=None, **kwargs):
    
        super().__init__(data, **kwargs)
        
        self.ID_COL = ID_COL
        self.LABEL_FREE = LABEL_FREE
        self.LABEL_GROUP_CAPTURE = LABEL_GROUP_CAPTURE
        self.PROTEASE = PROTEASE
        
    @property
    @cache 
    def unique_runs(self) -> set:
        if "Run" in self.columns: 
            return set(self["Run"].unique())
        else:
            return set()

    @property
    @cache
    def unique_ids(self) -> int:
        return self[self.ID_COL].nunique()


class LoadedData(pd.DataFrame):
    def __init__(self, loader):
        super().__init__(loader.data.copy())
        self.loader = loader

    def process(self, **kwargs) -> "ProcessedData":
        """
        Initiates processing. 
        Any kwargs passed here are forwarded 
        to the DataProcessor constructor.
        """
        from .processor import DataProcessor
        processor = DataProcessor(self.loader, **kwargs)
        return processor.process()