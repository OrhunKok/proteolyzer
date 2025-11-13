---
sidebar_label: cellenone
title: proteolyzer.cellenone.cellenone
---

## pd

## np

## os

## Literal

## Optional

## re

## MetaLogging

#### CELLEONE\_MAPPING

## CoordinatesMapping Objects

```python
class CoordinatesMapping(metaclass=MetaLogging)
```

#### \_\_init\_\_

```python
def __init__(root_dir: str,
             label_type: Optional[Literal["mTRAQ", "TMT"]] = None,
             plex: Optional[int] = None)
```

#### \_output\_file\_paths

```python
def _output_file_paths()
```

#### \_files\_parse

```python
def _files_parse(file_paths: dict)
```

#### \_data\_process

```python
def _data_process(parsed_data: dict)
```

#### \_metadata\_validate

```python
def _metadata_validate(metadata: pd.DataFrame)
```

#### map\_data

```python
def map_data() -> pd.DataFrame
```

#### \_stats\_process

```python
def _stats_process(parsed_stats: dict)
```

#### map\_stats

```python
def map_stats()
```

#### xls\_parse

```python
def xls_parse(file_paths: list)
```

#### log\_parse

```python
def log_parse(key: str, file_paths: list)
```

#### label\_well\_plex

```python
def label_well_plex(label_df: pd.DataFrame)
```

#### \_map\_coords

```python
def _map_coords(geo_df,
                map_df,
                coord_cols=["XPos", "YPos"],
                group_cols=["Target", "Field"])
```

