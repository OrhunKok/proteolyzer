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

## CHANNEL\_COL

## DROPLET\_COLS

## MERGE\_COLS

## NOZZLE\_WELL\_MAPPING

## PICKUP\_NOZZLE\_ID

## PICKUP\_NOZZLE\_XPOS\_OFFSET

## PRIMARY\_CHANNEL

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

#### \_classify\_file

```python
@staticmethod
def _classify_file(dirpath: str, filename: str) -> str | None
```

Which kind of output `filename` is, or None to ignore it.

The run directory name counts as well as the file name. Operators name
the logs inconsistently within one prep (&quot;labels_1&quot;, &quot;labeling_2&quot;,
&quot;L9&quot;, &quot;l4&quot;), while the ``.Run`` directory cellenONE creates keeps the
step name, so classifying on the file name alone silently files
labelling logs as dispense logs.

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

#### \_collapse\_label\_dispenses

```python
def _collapse_label_dispenses(label_df: pd.DataFrame) -> pd.DataFrame
```

One row per labelled position, carrying the droplets delivered to it.

A repeated labelling run dispenses the same channel onto a position
again, so the log holds several rows for it. Keeping them all
multiplies the merged metadata, so the latest dispense represents the
position and the droplet counts are summed. The raw ``Drops`` text
(&quot;58 drops&quot;) is replaced by two numbers:

Droplets
    Total droplets delivered to the position.
Dispenses
    How many dispense events delivered them.

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

#### \_merge\_imaging\_channels

```python
def _merge_imaging_channels(df: pd.DataFrame,
                            sample_name: str) -> pd.DataFrame
```

Collapse the per-imaging-channel rows of a cell onto one row.

cellenONE writes one geoprops row per (cell, channel): the
Transmission row carries the geometry, and each fluorescence channel
adds a row whose measurements are zero where nothing was detected.
Counting those as separate cells doubles the cell count and halves
every geometry average, so the extra channels become extra columns
(``Diameter.Green``) instead of extra rows.

#### \_cell\_key

```python
@staticmethod
def _cell_key(df: pd.DataFrame) -> list[str] | None
```

Columns identifying one cell, such that (cell, channel) is unique.

#### \_per\_channel\_columns

```python
@staticmethod
def _per_channel_columns(df: pd.DataFrame, key: list[str]) -> list[str]
```

Columns measured per channel, i.e. those differing within a cell.

Derived from the data rather than hardcoded, so an export with more or
fewer measurement columns still splits correctly.

#### log\_parse

```python
def log_parse(key: str, file_paths: list)
```

#### \_pickup\_plate

```python
@staticmethod
def _pickup_plate(file_path: str) -> int
```

Plate number for a pickup log, from its name or its run directory.

A pickup log is recognized by either, so the number has to be looked
for in both: `P7_..._Logfile.log` inside `pickup_7_....Run` carries it
only in the directory.

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

