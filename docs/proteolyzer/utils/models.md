---
sidebar_label: models
title: proteolyzer.utils.models
---

## pd

## List

## Union

## Set

## Path

## datetime

## logging

## re

## BaseModel

## Field

## field\_validator

## computed\_field

## Optional

## cached\_property

## Constants

## Data Objects

```python
class Data(BaseModel)
```

#### file\_path

#### load\_all\_columns

#### extra\_cols\_to\_load

#### \_validate\_extra\_cols\_to\_load

```python
@field_validator("extra_cols_to_load", mode="before")
@classmethod
def _validate_extra_cols_to_load(
        cls, value: Union[str, List[str], Set[str]]) -> set
```

#### \_validate\_path

```python
@field_validator("file_path", mode="after")
@classmethod
def _validate_path(cls, value: Path) -> Path
```

#### file\_name

```python
@computed_field
@cached_property
def file_name() -> str
```

#### file\_extension

```python
@computed_field
@cached_property
def file_extension() -> str
```

#### input\_type

```python
@computed_field
@cached_property
def input_type() -> str
```

#### cols\_subset

```python
@computed_field
@cached_property
def cols_subset() -> dict
```

#### cols\_rename\_mapping

```python
@computed_field
@cached_property
def cols_rename_mapping() -> dict
```

#### file\_stats

```python
@computed_field
@cached_property
def file_stats() -> dict
```

## ProcessedData Objects

```python
class ProcessedData(BaseModel)
```

#### data

#### ID\_COL

#### LABEL\_FREE

#### LABEL\_GROUP\_CAPTURE

#### PROTEASE

## Config Objects

```python
class Config()
```

#### arbitrary\_types\_allowed

Allow Pandas DataFrames

#### \_\_init\_\_

```python
def __init__(processor: "DataProcessor", **kwargs)
```

#### unique\_runs

```python
@computed_field
@cached_property
def unique_runs() -> set
```

#### unique\_ids

```python
@computed_field
@cached_property
def unique_ids() -> int
```

