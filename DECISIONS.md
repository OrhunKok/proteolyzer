# Decisions

Why this package is built the way it is. Most of these are already argued in the
docstring beside the code — this file is the index to them, so that a reader who
does not know the code exists can still find the reasoning. The ones marked
**measured** were settled with numbers rather than taste.

## Reading

**A format is a block on `core.formats.Config`, and nothing else knows its name.**
Detection walks whatever blocks are there, so adding an engine is a block and no
other change. It used to name DIANN and MaxQuant, which is why adding JMod and
FragPipe touched `models.py` at all.

**`LOAD_COLS` is per file, and a subset is an intersection.** A search writes
what its workflow produced: no IonQuant means no `Intensity`, no timsTOF means no
`Ion Mobility`, and JMod's parquet has seventeen of the columns its csv has. So
the loader intersects what was asked for with what the file has, and a consumer
checks a column arrived rather than assuming.

**`cols_to_load` replaces the subset; `extra_cols_to_load` adds to it.** The
configured subset is a floor, which is right for a pipeline and wrong for a
consumer that wants four columns of a large report. Precedence:
`load_all_columns`, then `cols_to_load`, then the configured subset plus extras.

**`rename=False` exists for consumers written against the engine's own names.**
This package normalizes every engine onto one schema, which is the right default
and the wrong one for a dashboard whose twenty-five pages say `Retention time`
and `Modified sequence`. Reversing the mapping afterwards is what it was doing.

**A file with no configured subset honours the asked-for columns.** It used to
read everything, so asking for two columns of an `allPeptides.txt` read all sixty
of them, and that file runs to several GB.

## Narrowing

**Categoricals are decided by measuring both representations.** *(measured)* Not
by a cardinality ratio, which is a poor proxy: on a real report the protein, gene
and name columns hold 0.3–0.4 distinct values per row and still halve in size,
while identifier columns approach 1.0, save nothing, and can come out *larger*.
`MIN_CATEGORICAL_SAVING = 0.2` separates the two groups, which measured 49%+ and
under 2%.

**A number is never made categorical**, however few values it holds: a
categorical of numbers no longer supports arithmetic, so a q-value column would
stop working downstream.

**Integers before floats, and both after the derived columns.** `RT.Width` is a
difference between two nearly equal retention times, so narrowing its inputs
first amplifies their rounding by the ratio between the time and the width — 6e-8
becomes 3e-5. A caller that derives its own columns afterwards passes
`narrow_floats=False`; `streamlit-DO-MS` does.

**Integrality is tested exactly.** `np.allclose` accepts a deviation of ~1000 at
the 1e8 magnitudes of an intensity column, which silently rounded quantitative
values.

**`round_large_floats` is off by default.** It discards real precision at the low
end of a quantitative column and saves nothing: the loss is below float32
resolution only above ~1.7e7, and the nullable integer it produces is wider than
the float32 it replaced.

**`Narrower` and `narrow(frame)` exist so the narrowing is usable without the
pipeline.** `process()` also renames, derives, labels and counts miscleavages,
and a consumer that does its own deriving wants none of that. Reaching the four
steps used to mean building a `Report` from a `Data` from a source that was never
read.

## cellenONE

**A geoprops row is a cell per imaging channel.** Counting rows doubles the cell
count and halves every geometry average. The transmission row is the cell and
every other channel joins onto it under its own name — `Diameter.Green` beside
`Diameter` — because a preparation that sorts on fluorescence sorts on the half
that filtering to transmission throws away.

**Which file is which step is worked out from the file.** Names are unreliable:
folders are numbered as often as named, and a directory that has been worked in
holds test runs, second attempts, and steps named after their reagent. The run's
own `Run ID:` line is read before its contents; the method summary next to it is
not, because a quench log names the label it is quenching.

**A parser starts the file at its start.** A run is read to work out which step
each file is and read again to parse it, and a consumed buffer parses as an empty
file.

## Packaging and release

**Optional subpackages are imported on first attribute access.** A bare core
install must not pay for matplotlib, seaborn or sqlalchemy, and
`tests/test_package_boundaries.py` runs a fresh interpreter to keep that true.

**Every tag gets a wheel.** A GitHub source archive carries no `.git`, so
setuptools-scm falls back and every tag installs as `0.0.0` — pip cannot then
tell two releases apart or compute an upgrade. `fallback_version` exists so an
archive *builds*; the wheel is what to pin.

**No LICENSE, deliberately.** This is public and released, and it still has no
license, which reads like an oversight and is not one: nothing on the account
gets one yet, because everything is early enough that the choice is not worth
making and would only have to be revisited. What it forecloses is publishing
outside the account — PyPI, or a dependant that is not private. Both consumers
are private and pin a wheel by URL, so neither is waiting on it. Revisit when
something is actually to be published, and not before.

**An invariant test is worth more than a case.** "Every engine block describes
only files it lists" found `msScans` in `LOAD_COLS` and missing from `FILES`,
which meant a `msScans.txt` was not recognised as MaxQuant output at all and came
back as one undelimited column. No example-based test was going to ask that.
