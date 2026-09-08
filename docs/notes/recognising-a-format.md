# Recognising a format

A file is identified in three passes, cheapest and most specific first, in
`core.models.Data.input_type`:

1. **By name.** Whatever a block lists in `FILES`, or matches with
   `FILE_PATTERNS`. A file called what its engine calls it is claimed without
   being opened.
2. **By whole shape.** `NARROW_SIGNATURES` — every declared column present, in a
   table about as narrow as declared.
3. **By an overlap of names.** `COLUMN_SIGNATURE` — any two of a set.

## Why the columns and not just the name

Because a file name is a convention people depart from. They rename what they
download, and one format has no default name to depart from in the first place:
Spectronaut's export is named by whoever ran the analysis, so
`GluC-30min.parquet` is as real a report as any. Before v0.11.0 that came back
`Unknown` — unrenamed, with a warning — although nothing about the file had
changed.

It is the same call the cellenONE reader had already made, for the same reason and
first: *names are unreliable*. Folders there are numbered as often as named, and a
directory that has been worked in holds test runs and second attempts.

## The four names that can be in no signature

Detection **refuses** a file two formats claim. That cuts both ways: a signature
reaching another engine's table would not mis-read one file, it would make *both*
unreadable. So before writing any signature, the real per-table column sets were
checked against each other, and four names are shared:

| name | shared by |
| --- | --- |
| `Charge` | MaxQuant and FragPipe |
| `Intensity` | MaxQuant and FragPipe |
| `PEP` | MaxQuant and DIA-NN |
| `rt` | JMod and DIA-NN's XIC export |

`Charge` and `Intensity` together would have made every FragPipe `psm.tsv`
ambiguous. That is the reason `test_every_real_table_is_claimed_by_exactly_one_format`
walks all twelve real tables of all five engines rather than trusting that nobody
reached for a shared name — it is the check that catches the fifth such name when
a sixth format arrives.

The column sets it checks against are written down in the test rather than derived
from the signatures, because a signature checked against itself checks nothing.
They came from the real DIA-NN report in `examples/`, from a real Spectronaut
export, and from the per-table lists `streamlit-DO-MS` keeps — which were
themselves recovered from the `LOAD_COLS` this package carried through v0.3.x.

Margins as measured: 3 to 15 columns matched against a threshold of 2.

## Why a narrow table needs its whole shape

DIA-NN's XIC export is `pr`, `feature`, `rt`, `value` — four words that belong to
nobody, one of which is JMod's. v0.11.0 recorded it as unsignable and left it
name-matched. That was giving up too early. What separates the two is not *which*
of those names is present but that **all** of them are and there is almost nothing
else:

```
xic columns:            ['feature', 'pr', 'rt', 'value']   (4)
shared with JMod:       ['rt']
JMod holds all four?    False
JMod width:             34
```

Either half separates this pair alone, which is the argument for requiring both.
Requiring the whole of a small schema is what works against another engine sharing
a name; the width is the second lock, against a frame carrying those four words
among two hundred others.

## Two limits, both deliberate

**A matching name still wins.** A FragPipe table renamed to `report.tsv` reads as
DIA-NN, because the name is asked first and answers. Preferring the columns would
mean opening every file to find out, and nothing here can tell a rename from a
report. The cost is a wrong rename mapping on a file somebody deliberately
misnamed.

**Two signatures matching is not an error.** A frame carrying two engines' columns
— somebody's DIA-NN report joined to a MaxQuant table for a figure, saved, read
back — warns and reads as `Unknown`. Two blocks claiming one *name* is this
package's config contradicting itself and does raise; a file that is nobody's is
what `Unknown` already means. v0.11.0 raised there by accident and v0.13.0 put it
back.

## What it costs

Nothing on the common path: the peek only happens when no block claims the name.
A parquet footer carries the schema and a delimited file gives up its header in
one line. `test_a_name_that_matches_is_still_taken_at_its_word` monkeypatches the
peek to *raise*, so a regression that started opening every file fails rather than
merely getting slower.

A peek that fails is answered with no columns rather than an exception — choosing
a reader is not the place to raise about a file, and the reader that follows will,
saying what it was actually trying to do.
