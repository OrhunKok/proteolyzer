# Spectronaut

Everything below was measured off real exports, not inferred from a specification.
Two of these facts shipped wrong first, which is the reason the note exists.

## The shape of a report

One row a precursor a run, long format. Columns are prefixed by the level they
belong to: `E.` the experiment, `R.` a run, `PG.` a protein group, `PEP.` a
peptide, `EG.` an elution group, `FG.` a fragment group — which is a precursor.

A report is configurable column by column, so what one lab's export holds is not
what another's does. That makes the intersection `cols_to_load` already takes
load-bearing rather than convenient: naming a column the analysis did not write
must not fail the read.

The two exports this was written from:

| | rows | runs | columns |
| --- | --- | --- | --- |
| tab separated, CRLF, 174 MB | 173,443 | 13 | 78 |
| parquet, 61 MB | 170,795 | 12 | 78 |

## There is no default output name

Whoever runs the analysis names the export. `<date>_<time>_<analysis>_Report` is
a shape an export often has and never has to — `GluC-30min.parquet` is a real
one. `FILE_PATTERNS` is kept as a shortcut for a report that does follow the
convention, but `COLUMN_SIGNATURE` is what actually identifies the format. See
[Recognising a format](recognising-a-format.md).

v0.6.0 leaned on the name alone, which is why v0.11.0 exists.

## The separator depends on the serialization

The tab-separated export writes `R.FileName`. **The parquet export writes
`R_FileName`**, and turns the space in `PG.Cscore (Run-Wise)` into an underscore
as well — a dot is a path separator in a nested parquet schema, so the export
spells it out of the way.

v0.8.0 added `.parquet` support on the stated assumption that the two agreed. They
do not, so all seventeen mapped columns missed and a report came back unrenamed:

```
mapping keys the real parquet has:              []
mapping keys present once '.' -> '_':  all 17 of them
```

Both spellings are carried now, derived from one list so they cannot drift.

## Columns written as text that are not text

Spectronaut is inconsistent about this **within one file**:

| column | stored as | should be |
| --- | --- | --- |
| `EG.Qvalue` | `'1.9967115436590949E-13'` | a float |
| `PG.Qvalue` beside it | a double | — |
| `PEP.IsProteotypic` | `'False'` | a bool |
| `EG.IsDecoy` beside it | a bool | — |
| `PEP.NrOfMissedCleavages` | `'1'` | an int |
| `EG.MinChannelQvalue`, `EG.MaxChannelQvalue` | scientific notation as text | floats |

A string q-value is not a lesser q-value; it raises `TypeError` the first time
anyone filters on it, which is the first thing anyone does with a report.

**Two columns are deliberately left as text.** `FG.XICDBID` is a database key
whose every value parses as a number — turning an identifier into an integer is
quiet damage, and it is the reason the retyped columns are a list and not a rule.
`EG.IsVerified` is every-value-`'NaN'` in the measured export, so there is nothing
there to say what it would be.

Of the export's 78 columns, exactly five are text that parses as numeric, and one
of those five is that key.

## `EG.InSourceFragmentationClass` is three states, not two and a gap

| value | rows | meaning |
| --- | --- | --- |
| `Likely Parent` | 1,117 | fragmented in the source |
| `Likely Child` | 1,146 | a fragment of one, with its parent named |
| `None` | 168,532 | neither |

The pairs are traceable:

```
Likely Parent  VVEAHVDQKNKVVTTPAFMCE    z=3  RT=26.979
Likely Child      AHVDQKNKVVTTPAFMCE    z=3  RT=26.982
    its parent: _VVEAHVDQKNKVVTTPAFMC[Carbamidomethyl (C)]E_.3
    co-elutes within 0.0034 min
```

The child is the parent minus its N-terminal `VVE` — a peptide losing residues in
the source and being identified as a shorter precursor, 0.2 seconds away. So
`None` there is a **category** and not an absence, and nulling it would erase a
real finding for 99% of the report. This is what
[Numbers and gaps](numbers-and-gaps.md) is about.

## Two more things not to assume

**`FG.Quantity` spans 2.54 to 400,000 in the one column**, so it is not on the
scale a DIA-NN area is. `round_large_floats` would throw away a fifth of a
precursor quantified at 2.54. It is off by default for every format, and there is
a test saying so about this one.

**`FG.PrecWindowNumber` is which of the method's isolation windows took the
precursor** — 60 of them on the measured export, agreeing with the same
instrument method's `--export-windows` design to within a Th. No other report read
here states it, so grouping by it recovers the window scheme with no design file
to hand. It is an integer, and a number is never made categorical, so nothing has
to keep it out of that.

## No precursor identifier

`EG.PrecursorId` is not in every export and was in neither measured one, so
`Precursor.Id` is built from `EG.ModifiedSequence` and `FG.Charge` — which is what
`streamlit-DO-MS` was doing by hand. It has to happen as the file is read, because
`DataProcessor` asks for that column in its constructor, before one of its own
steps runs.

Built only when renaming: those names are the core's vocabulary, and a caller
keeping the file's own has not asked for it.
