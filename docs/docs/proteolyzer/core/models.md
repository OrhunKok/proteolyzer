---
sidebar_label: models
title: proteolyzer.core.models
---

Typed models describing proteolyzer inputs and outputs.

:class:`Data` describes *where* data comes from and *what* should be read from
it; :class:`LoadedData` and :class:`ProcessedData` are thin ``DataFrame``
subclasses that carry the metadata needed by the next pipeline stage.

## datetime

## logging

## dataclass

## replace

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
def load() -> Report
```

Read the source into memory.

## Processing Objects

```python
@dataclass(frozen=True)
class Processing()
```

How a frame was processed, and what that revealed about it.

#### id\_col

Column identifying a precursor, used for the identification count.

#### label\_free

No labelling groups were found in the identifiers.

#### label\_group\_capture

Regex capturing a labelling group out of an identifier.

#### protease

Protease the missed-cleavage flag was computed for.

#### labels\_complete

False when channel information could not be derived in full.

#### rounded\_large\_floats

Large-magnitude float columns were rounded to integers.

## Report Objects

```python
@dataclass(frozen=True)
class Report()
```

A frame of proteomics data, with where it came from and what was done.

Composition rather than a ``DataFrame`` subclass: pandas returns plain
frames from most operations, so a subclass silently loses its metadata the
first time anyone slices it, and the pandas internals it has to hook into
are not a stable API. Use :attr:`frame` for anything pandas; the few
pass-throughs below are for interactive work.

#### frame

The data. Anything pandas goes through here.

#### source

What the frame was read from.

#### processing

Set once :meth:`process` has run.

#### is\_processed

```python
@property
def is_processed() -> bool
```

#### columns

```python
@property
def columns() -> pd.Index
```

#### \_\_len\_\_

```python
def __len__() -> int
```

#### \_\_getitem\_\_

```python
def __getitem__(key)
```

Column access, for interactive use. Returns a plain pandas object.

#### \_repr\_html\_

```python
def _repr_html_() -> str
```

#### process

```python
def process(**kwargs) -> Report
```

Normalize the frame: dtypes, derived columns, labelling information.

Keyword arguments are passed to :class:`~proteolyzer.core.processor
.DataProcessor`. Returns a new Report; this one is unchanged.

#### matrix

```python
def matrix(values: str, index: list[str], columns: list[str])
```

Pivot to a quantitative matrix. Returns a MatrixBuilder to chain on.

#### runs

```python
@property
def runs() -> set
```

The distinct runs present, empty if the frame has no Run column.

#### n\_identifications

```python
@property
def n_identifications() -> int
```

Distinct precursors, or 0 before processing.

#### with\_frame

```python
def with_frame(frame: pd.DataFrame) -> Report
```

The same report around a different frame, keeping its metadata.

