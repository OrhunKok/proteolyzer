---
sidebar_label: cellenone
title: proteolyzer.cellenone.cellenone
---

## os

## re

## Literal

## np

## pd

## Logged

## CELLEONE\_MAPPING

## DROPLET\_COLS

## MERGE\_COLS

## NOZZLE\_WELL\_MAPPING

## PICKUP\_NOZZLE\_ID

## PICKUP\_NOZZLE\_XPOS\_OFFSET

## TEMP\_STATS\_COLS

## CoordinatesMapping Objects

```python
class CoordinatesMapping(Logged)
```

#### \_\_init\_\_

```python
def __init__(root_dir: str,
             label_type: Literal["mTRAQ", "TMT"] | None = None,
             plex: int | None = None)
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
def _map_coords(
    geo_df,
    map_df,
    coord_cols=("XPos", "YPos"),
    group_cols=("Target", "Field")
) -> tuple[pd.Series, pd.Series]
```

Match each row of `map_df` to its nearest `plex` rows in `geo_df`.

Returns the matched geoprops indices and their distances, both as
Series of arrays indexed like `map_df` so they can be assigned straight
back onto it. With `plex` unset every candidate is kept, and the caller
resolves cells claimed more than once by distance.
