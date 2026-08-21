---
sidebar_label: processor
title: proteolyzer.core.processor
---

Higher-level processing pipelines and orchestrators.

This module composes loader, transformer and model utilities into
reusable processing functions and light-weight pipelines.

## warnings

## Callable

## cast

## np

## pd

## reference

## Config

## Logged

## Processing

## Report

## per\_distinct

#### CONFIG

## DataProcessor Objects

```python
class DataProcessor(Logged)
```

Processes raw data into a structured DataFrame.

#### \_\_slots\_\_

#### \_\_init\_\_

```python
def __init__(report: Report,
             id_col: str = "Precursor.Id",
             label_group_capture: str = r"\(((?:mTRAQ|SILAC|TMT)[^()]*)\)",
             protease: str = "Trypsin",
             round_large_floats: bool = False,
             narrow_floats: bool = True)
```

Initializes the DataProcessor.

Parameters
----------
report
    The report to process. Its frame is copied, so the input is
    untouched.
round_large_floats : bool, default False
    Round float columns whose median exceeds
    ``Config.COL_MEDIAN_THRESHOLD`` to integers, discarding their
    fractional part. Off by default because it loses real precision in
    the low range of quantitative columns.
narrow_floats : bool, default True
    Narrow float64 columns to float32. See
    :meth:`narrow_float_columns`; pass False to keep double precision.

#### process

```python
def process(verbose: bool = False) -> Report
```

Processes the data and returns a new :class:`Report`.

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

Drops columns holding the same value in every row.

Missing counts as a value, which reverses two cases that ``nunique``
alone gets backwards: it ignores NA, so a column that is entirely
empty counts 0 distinct values and survives, while one holding a
single value in a handful of rows counts 1 and is dropped -- losing
which rows those were.

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

#### narrow\_float\_columns

```python
def narrow_float_columns(df: pd.DataFrame) -> pd.DataFrame
```

Narrows float64 columns to float32 where the values allow it.

DIA-NN stores these columns as float32 in its own parquet output, so a
report read from parquet is already single precision while the same
report read from the TSV is not -- the identical data costs 44% more
memory depending on which file it came from. This makes the text path
match, at a worst-case relative error of 6e-8, float32&#x27;s epsilon,
measured across every float column of a real report.

Runs after :meth:`convert_float_columns_to_int`, so columns holding
whole numbers are already integers and keep their exact values;
``Ms1.Total.Signal.*`` reach 1e10, far beyond the 2**24 that float32
can represent exactly. It runs after :meth:`extra_info` too, so
derived columns are computed before anything is narrowed.

A column is left alone if any value falls outside float32&#x27;s normal
range, where narrowing would turn it into an infinity or a zero rather
than round it.

The bound is on each value, not on what is later computed from it:
subtracting two nearly equal narrowed values amplifies their rounding
by the ratio between them and their difference. Pass
``narrow_floats=False`` when doing that kind of arithmetic on the
frame directly.

#### narrow\_integer\_columns

```python
def narrow_integer_columns(df: pd.DataFrame) -> pd.DataFrame
```

Narrows integer columns to the smallest dtype that holds them.

Unlike the float narrowing this is exact: an integer either fits or it
does not, and the width is chosen per column. A charge state, a run
index or a peptide length arrives as int64 and needs one byte, while
``Ms1.Total.Signal`` reaches 1e10 and keeps all eight.

Runs after :meth:`extra_info` so the integer columns derived there are
narrowed too. Nullable dtypes stay nullable.

The values are unchanged, but the headroom is not: adding a scalar
that takes a narrow column past its range wraps around rather than
raising, as it did before for the columns
:meth:`convert_float_columns_to_int` narrows. Reductions such as
``sum`` accumulate in int64 and are unaffected.

#### \_fits\_in\_float32

```python
@staticmethod
def _fits_in_float32(column: pd.Series) -> bool
```

Whether every value of `column` is within float32&#x27;s normal range.

#### convert\_columns\_to\_categorical

```python
def convert_columns_to_categorical(df: pd.DataFrame) -> pd.DataFrame
```

Converts columns to categorical where that actually saves memory.

Decided by measuring both representations rather than by a cardinality
ratio, which is a poor proxy for the saving: on a real report the
protein, gene and name columns hold 0.3-0.4 distinct values per row
and still halve in size, while identifier columns approach 1.0, save
nothing, and can come out *larger* than the strings they replace. A
ratio tight enough to exclude the second group excludes most of the
first as well.

Numeric columns are never converted, however few distinct values they
hold: a categorical of numbers no longer supports arithmetic, so
q-value or intensity columns would stop working downstream.

#### \_as\_categorical

```python
@staticmethod
def _as_categorical(column: pd.Series) -> tuple[float, pd.Series]
```

`column` as a categorical, and the fraction of memory that saves.

Measured by building it rather than estimated: estimating means
guessing the width of the codes and the size of the dictionary, and
the factorize it costs is what counting distinct values would cost
anyway. The conversion is handed back so a caller that decides to keep
it does not build it a second time. The saving is negative when the
categorical is the larger of the two.

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

Every column derived here is a function of the precursor identifier alone,
so the regex extraction and the pivots that follow it run over the
*distinct* identifiers and the result is gathered back out. A report holds
one row per identifier per run per channel, so on a 40-run experiment the
regex extraction -- which pandas applies element by element in Python, on
categorical columns too -- does a fortieth of the work.

#### \_\_slots\_\_

#### \_\_init\_\_

```python
def __init__(processed_data: DataProcessor)
```

Initializes the LabelGenerator.

#### \_distinct\_ids

```python
def _distinct_ids() -> tuple[pd.Series, np.ndarray]
```

The distinct identifiers, and which one each row holds.

Keyed by position rather than by the identifier itself: the frames
derived from these get a plain integer index, so the pivots below sort
a few thousand integers instead of a few hundred thousand long
strings, and attaching the result is a positional gather rather than a
join on those strings.

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

Attaches the per-identifier label information to the data.

The offsets are derived last: they need the counts to average over,
and they convert a column of `sorted_matches` in place.

#### \_generate\_run\_channels

```python
def _generate_run_channels(df: pd.DataFrame) -> pd.DataFrame
```

Generates run channel information.

