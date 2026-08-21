---
sidebar_label: models
title: proteolyzer.utils.models
---

Typed models describing proteolyzer inputs and outputs.

:class:`Data` describes *where* data comes from and *what* should be read from
it; :class:`LoadedData` and :class:`ProcessedData` are thin ``DataFrame``
subclasses that carry the metadata needed by the next pipeline stage.

## datetime

## logging

## cached\_property

## Path

## IO

## pd

## BaseModel

## ConfigDict

## Field

## computed\_field

## field\_validator

## Config

#### CONFIG

#### logger

#### SourceType

## Data Objects

```python
class Data(BaseModel)
```

A description of a single input file (or file-like object).

#### model\_config

#### source

#### load\_all\_columns

#### extra\_cols\_to\_load

#### INPUT\_TYPE

#### is\_path

```python
@property
def is_path() -> bool
```

#### is\_file\_like

```python
@property
def is_file_like() -> bool
```

#### \_path

```python
@property
def _path() -> Path
```

The source as a Path. Only meaningful when :attr:`is_path`.

#### \_validate\_source

```python
@field_validator("source", mode="after")
@classmethod
def _validate_source(cls, value: SourceType) -> SourceType
```

#### \_validate\_extra\_cols\_to\_load

```python
@field_validator("extra_cols_to_load", mode="before")
@classmethod
def _validate_extra_cols_to_load(cls, value: object) -> set[str] | None
```

#### file\_name

```python
@computed_field
@cached_property
def file_name() -> str
```

Stem of the source file, used to recognize known input formats.

#### file\_extension

```python
@computed_field
@cached_property
def file_extension() -> str
```

#### file\_stats

```python
@computed_field
@cached_property
def file_stats() -> dict | None
```

#### input\_type

```python
@computed_field
@cached_property
def input_type() -> str
```

The search engine that produced this file: DIANN, MaxQuant or Unknown.

#### cols\_subset

```python
@computed_field
@cached_property
def cols_subset() -> set[str] | None
```

Columns to read, or ``None`` to read everything.

#### cols\_rename\_mapping

```python
@computed_field
@cached_property
def cols_rename_mapping() -> dict
```

#### load

```python
def load() -> "LoadedData"
```

Read the source into memory.

## ProcessedData Objects

```python
class ProcessedData(pd.DataFrame)
```

A specialized DataFrame to hold processed data and its metadata.
Inherits all pandas DataFrame methods and attributes.

#### \_metadata

#### \_constructor

```python
@property
def _constructor()
```

#### \_\_init\_\_

```python
def __init__(data=None,
             ID_COL=None,
             LABEL_FREE=None,
             LABELS_COMPLETE=None,
             LABEL_GROUP_CAPTURE=None,
             PROTEASE=None,
             **kwargs)
```

#### unique\_runs

```python
@property
def unique_runs() -> set
```

#### unique\_ids

```python
@property
def unique_ids() -> int
```

## LoadedData Objects

```python
class LoadedData(pd.DataFrame)
```

Raw data as read from disk, plus the loader that produced it.

#### \_metadata

#### \_\_init\_\_

```python
def __init__(loader)
```

#### process

```python
def process(**kwargs) -> ProcessedData
```

Initiates processing.
Any kwargs passed here are forwarded
to the DataProcessor constructor.
