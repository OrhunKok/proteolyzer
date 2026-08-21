---
sidebar_label: loader
title: proteolyzer.utils.loader
---

Data loading helpers for common proteolyzer inputs.

This module implements robust, documented helpers to read CSV/TSV, Excel,
and other domain-specific export formats. Helpers validate required
columns and return pandas DataFrame objects.

## csv

## Callable

## partial

## Path

## chardet

## pd

## pq

## Config

## Logged

## Data

#### CONFIG

## DataLoader Objects

```python
class DataLoader(Logged)
```

Loads data from various file formats with robust handling.

#### \_\_slots\_\_

#### \_\_init\_\_

```python
def __init__(file: Data, verbose: bool = False)
```

Initializes DataLoader.

Parameters
----------
file : Data
    The source descriptor. `file.source` may be a Path, string path,
    or file-like object. It is held rather than copied field by field,
    so the two cannot drift apart.
verbose : bool, default False
    If True, run optional diagnostic steps after loading (currently:
    log the in-memory size of the loaded data). Off by default because
    memory_usage(deep=True) walks every cell of every object column,
    which is expensive on large frames.

#### source

```python
@property
def source()
```

#### INPUT\_TYPE

```python
@property
def INPUT_TYPE() -> str
```

#### cols\_rename\_mapping

```python
@property
def cols_rename_mapping() -> dict
```

#### cols\_subset

```python
@property
def cols_subset()
```

#### file\_name

```python
@property
def file_name() -> str
```

#### file\_extension

```python
@property
def file_extension() -> str
```

Lower-cased suffix; taken from Data so named streams dispatch too.

#### is\_path

```python
@property
def is_path() -> bool
```

#### \_auto\_load

```python
def _auto_load() -> pd.DataFrame
```

Automatically dispatch to the correct loader based on file extension.

#### \_rewind

```python
def _rewind() -> None
```

Return a file-like source to the start so it can be read again.

Peeking at the header and then reading the body means the source is read
twice; paths reopen transparently but streams have to be rewound.

#### \_rename\_cols

```python
def _rename_cols(df: pd.DataFrame | None = None) -> pd.DataFrame
```

#### \_cols\_to\_load

```python
def _cols_to_load(all_cols) -> list
```

Intersect the file&#x27;s columns with the configured subset.

Keeps the order the columns appear in the file, so repeated loads of the
same file always produce the same column order.

#### \_load\_csv

```python
def _load_csv(delimiter: str | None = None) -> pd.DataFrame
```

#### \_get\_delimiter

```python
def _get_delimiter(default_delimiter="\t",
                   sample_size=524288,
                   sample_percent=0.01) -> str
```

Detect delimiter for CSV-like files.

#### \_load\_excel

```python
def _load_excel() -> pd.DataFrame
```

#### \_load\_parquet

```python
def _load_parquet() -> pd.DataFrame
```

#### \_load\_txt

```python
def _load_txt() -> pd.DataFrame
```

Load a text file. Detect encoding and MIME type if Path.

#### \_detect\_encoding

```python
def _detect_encoding(sample_size: int = 100_000) -> str
```

Detect the character encoding of a file path.
