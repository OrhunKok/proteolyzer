# Numbers and gaps

Two questions a reader has to answer the same way whichever file it is looking at:
which cells are empty, and what a number is.

## What counts as missing

One definition, used by both readers, so that a report read as text and the same
report read as parquet agree about which cells are empty. Everything in it is a
machine artefact:

```
""  NaN  nan  -NaN  -nan  <NA>  #N/A  #N/A N/A  #NA  1.#IND  -1.#IND  1.#QNAN  -1.#QNAN
```

`NaN` is `str(float("nan"))`, `<NA>` is pandas' own, `1.#IND` and `#N/A` are what
a C runtime and a spreadsheet write. None is a word in any language, so nothing is
lost reading them as the gaps they are.

## Why no words

`NA`, `N/A`, `None`, `null` and `NULL` are all nulled by `read_csv` unless it is
told otherwise, and this package took that default. It should not have:
Spectronaut's `EG.InSourceFragmentationClass` is a three-state classification
whose third state is the word `None`, in 168,532 rows of 170,795 —
see [Spectronaut](spectronaut.md#eginsourcefragmentationclass-is-three-states-not-two-and-a-gap).
Taking the default erased that state into "not recorded" for 99% of the report,
and the same file read as parquet kept it, so the two serializations disagreed
about **data** rather than about dtype.

That a column *may* use one of these words for an absence is true. That it may use
one as a value is also true, and only one of those two mistakes is silent: the one
this makes shows up in the dtype. A numeric column that spells its holes `NA`
rather than leaving them empty comes back as text, which a caller notices. There
is a test that a numeric column with ordinary empty holes still reads numeric,
which is the obvious way to get this wrong.

The four other engines read here all leave their holes empty, so nothing changes
for them.

## Where the parquet path needed help

A delimited reader turns these into gaps itself, being handed the set. Parquet
stores the string it was given, so a pass afterwards is what makes the two agree —
and the parquet one lied the worse way round, because **`pd.isna` says `False` of
the string `"NaN"`**, so a gap read as data all the way to whatever plotted it.

## Columns written as text

Which columns a format writes as text is a fact about the format, so the format
block says. It is a **list and not a rule** because the rule gets it wrong:
`FG.XICDBID` is a database key whose every value parses as a number.

Conversion is **all or nothing** per column. `errors="coerce"` would hand back a
numeric column shorter by however many values it could not read, with nobody told;
a column where any value is not a number is left exactly as it arrived and logged
with why. A column an export already stored properly is not touched, so a lab
whose export got it right is not round-tripped for nothing.

Applied before the rename, on the file's own names, so it reaches a caller reading
with `rename=False` — a q-value that cannot be compared to a float is no more use
under one name than another.

## pandas loses digits

This one reaches further than the columns being retyped:

```
python float()            0.0007162974636774825
read_csv engine=pyarrow   0.0007162974636774825   ==
pyarrow cast              0.0007162974636774825   ==
read_csv engine=c         0.0007162974636774      !=
pd.to_numeric             0.0007162974636774      !=
```

pandas' float parser stops at about sixteen significant digits. A relative error
of 1.2e-13 is four orders past float64's epsilon, and it is truncation rather than
rounding.

`core.loader` says the stock parser is a fallback for when memory is short, giving
identical values, and a test says the fast parser must not change the frame.
Neither was true for a decimal that long — the fixture simply had no value long
enough to notice. `float_precision="round_trip"` makes them true, at ~3x on the
parse, on a path already chosen for fitting in memory rather than for being quick.
Both paths are pinned by tests that read one file each way and compare.
