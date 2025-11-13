---
sidebar_label: loader
title: proteolyzer.utils.loader
---

Data loading helpers for common proteolyzer inputs.

This module implements robust, documented helpers to read CSV/TSV, Excel,
and other domain-specific export formats. Helpers validate required
columns and return pandas DataFrame objects.

## pd

## partial

## csv

## magic

## chardet

## pq

## Data

## MetaLogging

## DataLoader Objects

```python
class DataLoader(metaclass=MetaLogging)
```

Loads data from various file formats.

#### \_\_slots\_\_

#### \_\_init\_\_

```python
def __init__(file: Data)
```

Initializes the DataLoader.

#### \_auto\_load

```python
def _auto_load() -> pd.DataFrame
```

Automatically loads data based on file extension.

#### \_rename\_cols

```python
def _rename_cols() -> pd.DataFrame
```

Renames columns based on the mapping.

#### \_get\_delimiter

```python
def _get_delimiter(default_delimiter: str = "\t",
                   encoding: str = "utf-8",
                   min_sample_size: int = 524288,
                   sample_percent: float = 0.01) -> str
```

Detects the delimiter of a CSV-like file.

#### \_cols\_to\_load

```python
def _cols_to_load(all_cols: set) -> list
```

Determines columns to load based on settings.

#### \_load\_csv

```python
def _load_csv(delimiter: str = None) -> pd.DataFrame
```

Loads a CSV or TSV file.

#### \_load\_excel

```python
def _load_excel() -> pd.DataFrame
```

Loads an Excel file.

#### \_load\_parquet

```python
def _load_parquet() -> pd.DataFrame
```

Loads a Parquet file.

#### \_load\_txt

```python
def _load_txt() -> pd.DataFrame
```

Loads a plaintext file with detected MIME type and encoding.

#### \_detect\_file\_type\_and\_encoding

```python
def _detect_file_type_and_encoding()
```

Strictly detects MIME type and text encoding for the file. Raises if undetectable.

