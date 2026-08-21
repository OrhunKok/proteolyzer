---
sidebar_label: loader
title: proteolyzer.core.loader
---

Data loading helpers for common proteolyzer inputs.

This module implements robust, documented helpers to read CSV/TSV, Excel,
and other domain-specific export formats. Helpers validate required
columns and return pandas DataFrame objects.

## csv

## os

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

#### BULK\_READ\_EXPANSION

Peak memory the fast CSV parser needs, as a multiple of the file&#x27;s size on
disk. Measured at ~16x on a wide, float-heavy report; rounded up, since
overshooting only costs time and undershooting costs the whole read.

#### ASSUMED\_AVAILABLE\_MEMORY

Assumed free memory where the platform will not say, chosen to be smaller
than any machine that would be running a search in the first place.

#### \_available\_memory

```python
def _available_memory() -> int
```

Free physical memory in bytes, or :data:`ASSUMED_AVAILABLE_MEMORY`.

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

#### \_read\_delimited

```python
def _read_delimited(delimiter: str, cols_to_load: list) -> pd.DataFrame
```

Read the body, preferring pyarrow&#x27;s multithreaded CSV reader.

A search report is the largest thing proteolyzer reads, and on a
million-row report pyarrow is ~4x faster than the stock parser whole,
~11x when reading a subset of the columns, for identical dtypes and
values. It is stricter, though -- it rejects ragged rows the stock
parser pads -- so that one stays as a fallback.

It also builds an Arrow table and the frame at the same time, so it
needs several times the memory the frame ends up taking; past
:data:`BULK_READ_EXPANSION` it is not worth the risk of running out.

The two parsers agree on column order only because
:meth:`_cols_to_load` hands over the columns in the order the file has
them: pyarrow returns them in the order asked for, the stock parser
always in file order.

#### \_fast\_read\_fits

```python
def _fast_read_fits() -> bool
```

Whether the fast parser&#x27;s peak memory is affordable for this source.

A stream is already in memory, so there is nothing to weigh. For a
path, the estimate is the file&#x27;s size on disk: text expands as it is
parsed, and the Arrow table and the frame coexist, which measured
~16x the file size in peak RSS on a wide, float-heavy report.

Free physical memory is what is compared against, which underestimates
what is usable (page cache is reclaimable) and, inside a container,
reports the host&#x27;s rather than the cgroup&#x27;s. Both err towards the
parser that needs less memory, which costs time and nothing else.

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

