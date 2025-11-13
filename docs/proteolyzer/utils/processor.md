---
sidebar_label: processor
title: proteolyzer.utils.processor
---

## pd

## np

## warnings

## MetaLogging

## ProcessedData

## Constants

## DataLoader

## DataProcessor Objects

```python
class DataProcessor(metaclass=MetaLogging)
```

Processes raw data into a structured DataFrame.

#### \_\_slots\_\_

#### \_\_init\_\_

```python
def __init__(data_loader: DataLoader,
             ID_COL: str = "Precursor.Id",
             LABEL_GROUP_CAPTURE: str = r"\((?!.*\bUniMod\b)(.*?)\)",
             PROTEASE: str = "Trypsin")
```

Initializes the DataProcessor.

#### process

```python
def process() -> ProcessedData
```

Processes the data and returns a ProcessedData object.

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

Converts eligible float columns to integer type.

#### convert\_columns\_to\_categorical

```python
def convert_columns_to_categorical(df: pd.DataFrame) -> pd.DataFrame
```

Converts eligible columns to categorical type.

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

Calculates miscleavages and adds a boolean column.

## \_LabelGenerator Objects

```python
class _LabelGenerator(metaclass=MetaLogging)
```

Generates label information for DIA-NN data.

#### \_\_slots\_\_

#### \_\_init\_\_

```python
def __init__(processed_data: ProcessedData)
```

Initializes the LabelGenerator.

#### label\_checked\_data

```python
@property
def label_checked_data()
```

Returns the label-checked data.

#### \_validate\_matrix\_shape

```python
def _validate_matrix_shape(matrix: pd.DataFrame) -> pd.DataFrame
```

Validates the shape of the label matrix.

#### \_label\_matrix

```python
def _label_matrix(sorted_matches: pd.DataFrame) -> pd.DataFrame
```

Generates the label matrix.

#### \_label\_counts

```python
def _label_counts(sorted_matches: pd.DataFrame) -> pd.DataFrame
```

Generates the label counts matrix.

#### \_label\_offset

```python
def _label_offset(sorted_matches: pd.DataFrame,
                  label_count: pd.DataFrame) -> pd.DataFrame
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

