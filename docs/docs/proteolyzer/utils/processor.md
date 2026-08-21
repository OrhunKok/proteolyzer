---
sidebar_label: processor
title: proteolyzer.utils.processor
---

Higher-level processing pipelines and orchestrators.

This module composes loader, transformer and model utilities into
reusable processing functions and light-weight pipelines.

## warnings

## Callable

## cast

## np

## pd

## Config

## DataLoader

## Logged

## ProcessedData

#### CONFIG

#### \_protease\_rules

```python
def _protease_rules(protease: str)
```

Look up the cleavage rules for `protease`, with a helpful error if absent.

## DataProcessor Objects

```python
class DataProcessor(Logged)
```

Processes raw data into a structured DataFrame.

#### \_\_slots\_\_

#### \_\_init\_\_

```python
def __init__(data_loader: DataLoader,
             ID_COL: str = "Precursor.Id",
             LABEL_GROUP_CAPTURE: str = r"\(((?:mTRAQ|SILAC|TMT)[^()]*)\)",
             PROTEASE: str = "Trypsin",
             ROUND_LARGE_FLOATS: bool = False)
```

Initializes the DataProcessor.

Parameters
----------
ROUND_LARGE_FLOATS : bool, default False
    Round float columns whose median exceeds
    ``Config.COL_MEDIAN_THRESHOLD`` to integers, discarding their
    fractional part. Off by default because it loses real precision in
    the low range of quantitative columns.

#### process

```python
def process(verbose: bool = False) -> ProcessedData
```

Processes the data and returns a ProcessedData object.

Parameters
----------
verbose : bool, default False
    If True, run optional diagnostic steps after processing (currently:
    log the in-memory size of the processed data).

#### \_check\_labelfree

```python
def _check_labelfree() -> None
```

Checks if the data is label-free.

#### drop\_identical\_cols

```python
def drop_identical_cols(df: pd.DataFrame) -> pd.DataFrame
```

Drops columns with identical values.

#### convert\_float\_columns\_to\_int

```python
def convert_float_columns_to_int(df: pd.DataFrame) -> pd.DataFrame
```

Narrows float columns that already hold whole numbers to integers.

Integrality is tested exactly. ``np.allclose`` is not usable here: its
default relative tolerance of 1e-5 accepts a deviation of ~1000 at the
1e8 magnitudes typical of intensity columns, which silently rounded
quantitative values.

Columns that are *not* integral are left alone unless
``ROUND_LARGE_FLOATS`` is set; see :meth:`_round_large_floats`.

#### \_to\_narrow\_int

```python
@staticmethod
def _to_narrow_int(column: pd.Series, values: np.ndarray, complete: bool)
```

Return `column` as the smallest integer dtype that holds it exactly.

A plain numpy integer dtype is used when nothing is missing; otherwise a
nullable one, which needs an extra byte per value for its mask. Picking
the width matters: unconditionally using Int64 makes the float32 columns
DIA-NN writes *larger*, defeating the point of narrowing them.

#### \_round\_large\_floats

```python
def _round_large_floats(df: pd.DataFrame) -> pd.DataFrame
```

Optionally round large-magnitude float columns to integers.

Opt-in (``ROUND_LARGE_FLOATS``) because it discards real precision at
the low end of a quantitative column while saving nothing: the loss is
below float32 resolution only above ~1.7e7, and the resulting nullable
integer dtype is wider than the float32 it replaces.

#### convert\_columns\_to\_categorical

```python
def convert_columns_to_categorical(df: pd.DataFrame) -> pd.DataFrame
```

Converts eligible low-cardinality columns to categorical type.

Numeric columns are never converted, however few distinct values they
hold: a categorical of numbers no longer supports arithmetic, so
q-value or intensity columns would stop working downstream.

#### rename\_columns

```python
def rename_columns(df: pd.DataFrame) -> pd.DataFrame
```

Renames columns based on given alias mapping.

#### extra\_info

```python
def extra_info(df: pd.DataFrame) -> pd.DataFrame
```

Adds extra information columns.

#### miscleavages

```python
def miscleavages(df: pd.DataFrame,
                 seq_col: str = "Stripped.Sequence",
                 protease: str = "Trypsin") -> pd.DataFrame
```

Flags peptides whose residue counts are inconsistent with full cleavage.

## \_LabelGenerator Objects

```python
class _LabelGenerator(Logged)
```

Generates label information for DIA-NN data.

#### \_\_slots\_\_

#### \_\_init\_\_

```python
def __init__(processed_data: DataProcessor)
```

Initializes the LabelGenerator.

#### \_validate\_matrix\_shape

```python
def _validate_matrix_shape(matrix: pd.DataFrame | None,
                           categorical: bool = True) -> pd.DataFrame | None
```

Validates the shape of the label matrix.

`categorical` is off for the count matrix: it holds numbers, and a
categorical of numbers no longer supports arithmetic.

#### \_label\_matrix

```python
def _label_matrix(sorted_matches: pd.DataFrame) -> pd.DataFrame | None
```

Generates the label matrix.

#### \_label\_counts

```python
def _label_counts(sorted_matches: pd.DataFrame) -> pd.DataFrame | None
```

Generates the label counts matrix.

#### \_label\_offset

```python
def _label_offset(sorted_matches: pd.DataFrame,
                  label_count: pd.DataFrame | None) -> pd.DataFrame | None
```

Generates the label offset matrix.

#### \_generate\_sorted\_matches

```python
def _generate_sorted_matches() -> pd.DataFrame
```

Generates sorted matches DataFrame.

#### \_add\_label\_info

```python
def _add_label_info(df: pd.DataFrame,
                    sorted_matches: pd.DataFrame) -> pd.DataFrame
```

Adds label information to the DataFrame.

#### \_generate\_run\_channels

```python
def _generate_run_channels(df: pd.DataFrame) -> pd.DataFrame
```

Generates run channel information.
