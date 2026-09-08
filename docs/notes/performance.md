# Performance

Measured on the 170,795-row × 78-column Spectronaut parquet export, on a busy
machine — so every comparison below is **interleaved and best-of-N in one
process**, which is the only way to get a usable number at load average 13.

## Where a read's time goes

| step | time |
| --- | --- |
| `_load_parquet` — the file itself | 5.3 s for all 78 columns |
| `_retype_cols` | 1.19 s |
| `_close_text_gaps` | 0.21 s |
| `_rename_cols` | 0.007 s |

The whole read is ~1.3 s when a caller names the columns it wants, which is what
`cols_to_load` is for. Reading the file dominates; to go faster from here, read
fewer columns.

## Parsing a number: Arrow, not pandas

`_retype_cols` cost as much as reading the file, for five columns. Broken down
over one column of 170,795 values, the answer was not where it looked:

| step | time |
| --- | --- |
| `astype("string")` | 0.000 s |
| `.str.strip()` | 0.003 s |
| `.isin(...)` | 0.003 s |
| **`pd.to_numeric`** | **0.153 s** |
| **pyarrow cast** | **0.007 s** |

All the string handling was free; the parse was everything. Casting the Arrow
buffer the values already sit in is **23x** quicker — ~6x on the step and 1.34x on
a whole parquet read.

It is also the more accurate of the two, which turned out to matter more than the
speed: see [Numbers and gaps](numbers-and-gaps.md#pandas-loses-digits). `to_numeric`
remains the fallback for anything Arrow will not take.

## One optimisation measured and thrown away

Skipping `_close_text_gaps` for readers handed the gap set, which have already
done the work. The reasoning was right — that pass finds nothing after a text read
— and the measurement said:

```
text read, gap pass still run : 0.423s
text read, gap pass skipped   : 0.442s
```

No saving, because an `isin` over an Arrow string column is 3 ms; it is ~1.7% of a
read, below the noise floor. And it would have coupled the loader to which readers
get `na_values`, where drifting out of step brings back the exact bug the pass
exists to prevent. Not worth it. The comment on that line records the result so
nobody re-derives it.

## Two false alarms, for the next person who benchmarks here

**The suite went from 7 s to 35 s** mid-change and looked like a regression. It was
the machine: load average 11.9 on 10 cores. Profiling put the time in the parquet
read and in steps that predated the change.

**A sequential A/B showed 3.30 s against 1.14 s**, which was convincing and wrong:
the two arms ran one after the other on a warming page cache and a dropping load.
Interleaved, the same comparison came out 1.74 s against 2.00 s — i.e. no
difference. Never measure two arms in sequence on this box.
