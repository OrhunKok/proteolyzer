# Memory and dtypes

A search report is the largest thing proteolyzer handles, so `process()` works
to make one cheaper to hold. On real reports it comes down by a third to a
half:

| processed frame | before | after |
|---|---|---|
| `report.parquet` (10k × 67) | 4.93 MB | **3.74 MB** |
| `report.tsv` (the same data) | 6.59 MB | **3.75 MB** |
| 500k × 59 TSV | 234.98 MB | **120.75 MB** |

## Where it goes

```python
processed.memory()
```

| Column | Dtype | Bytes | Share |
|---|---|---|---|
| Precursor.Id | str | 334133 | 8.2% |
| Modified.Sequence | str | 324133 | 8.0% |
| Protein.Ids | category | 130726 | 3.2% |

Counted deeply, so the strings behind a text column are measured rather than
the pointers to them. This is the first place to look when a report is larger
than expected.

## Integers

Float columns that already hold whole numbers become the smallest exact integer
dtype, and integer columns narrow to the width they need. Both are exact: a
value either fits or it does not.

A charge state, a run index or a peptide length arrives as `int64` and needs one
byte. `Ms1.Total.Signal` reaches 1e10 and keeps all eight.

!!! warning "Narrow integers have less headroom"

    The values do not change, but the room above them does. Adding a scalar
    that takes a narrow column past its range wraps around rather than raising.
    Reductions such as `sum` accumulate in `int64` and are unaffected.

Columns with a real fractional part keep it. To also round high-magnitude
columns to integers — discarding their fractional part — ask for it:

```python
report.process(round_large_floats=True)
```

That is off by default because it loses real precision at the low end of a
quantitative column.

## Single precision

Remaining `float64` columns narrow to `float32`. That is what DIA-NN stores in
its own parquet output, so a report read from the TSV used to cost 44% more
than the same report read from parquet; now they cost the same.

The worst relative error measured across every float column of a real report is
**6e-8**, float32's epsilon. A column is left alone if any of its values falls
outside float32's normal range, where narrowing would substitute an infinity or
a zero rather than round it.

```python
report.process(narrow_floats=False)   # keep double precision
```

!!! warning "Differences of nearly equal values"

    The 6e-8 bound is on each value, not on what you compute from it.
    Subtracting two nearly equal narrowed values amplifies their rounding by
    the ratio between them and their difference. `RT.Width` is exactly that
    kind of quantity, so it is derived *before* anything is narrowed — but if
    you do that sort of arithmetic on the frame yourself, pass
    `narrow_floats=False`.

## Categoricals

Text columns become categorical when measuring both representations says that
saves at least 20% of the column.

Cardinality is a poor proxy for the saving. On a real report the protein, gene
and name columns hold 0.3–0.4 distinct values per row and still halve, while
identifier columns approach 1.0, save nothing, and can come out *larger* than
the strings they replace. Any ratio tight enough to exclude the second group
excludes most of the first as well.

Numeric columns are never converted, however few distinct values they hold: a
categorical of numbers no longer supports arithmetic, so a q-value or intensity
column would stop working downstream.

## Reading

The pyarrow parser is several times faster than the stock one — on a
million-row report, 4× reading it whole and 11× reading a subset of the columns
— for identical dtypes and values.

It also builds an Arrow table and the frame at the same time, needing roughly
16× the file's size on disk in peak memory, which is 1.8× the stock parser's
peak. So it is used only when the file is small enough for that against the
memory available; past that the stock parser reads it with a bounded footprint.
Both produce identical frames, so the choice trades only time.

Two differences are worth knowing, since both fall back to the stock parser:

- pyarrow rejects ragged rows that the stock parser pads with NA.
- A column pyarrow cannot decode as UTF-8 comes back as raw bytes rather than
  raising. Falling back means you get the stock parser's `UnicodeDecodeError`
  instead of a column of `bytes`.

## Deriving repeated values once

pandas applies string operations element by element in Python, even on
categorical columns, and a report repeats every identifier, peptide and protein
group once per run per channel.
[`per_distinct`][proteolyzer.core.operations.per_distinct] applies an
element-wise transformation to the distinct values and gathers the results back
out, which is the same answer for a fraction of the work — an order of
magnitude on a multi-run report.

```python
from proteolyzer.core.operations import per_distinct

razor = per_distinct(lambda x: x.str.split(";").str[0])(frame["Protein.Group"])
```

The transformation must be element-wise: the result for a row may depend on
that row's value and nothing else.
