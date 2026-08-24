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

**Which columns to keep is the caller's, and this package keeps no list.** A
format block says what the file *is* -- which names it goes by, how its columns
map onto the canonical schema, which must stay numeric. It used to say which
columns were worth reading as well, 260 of them across 29 tables, contributed
from the dashboard so that four repositories would not each carry their own.

That was the wrong seam. Which columns matter is a fact about the project doing
the reading, not about the file: a dashboard plots the m/z and the injection time
a pipeline never looks at, and a pipeline wants quantities no panel shows. The
tell was that the one consumer the lists came from had already overridden two of
them, with a comment saying the core's subset "is a pipeline's rather than a
dashboard's". A shared list that every sharer overrides is a shared list in name.

So the core reads the file whole and `Data.cols_to_load` is where a project says
otherwise. A reader that silently drops a column is worse than a wide frame: the
column is gone with nothing said, and the caller who needed it finds out from a
`KeyError` three functions later.

**A subset is an intersection.** A search writes what its workflow produced: no
IonQuant means no `Intensity`, no timsTOF means no `Ion Mobility`, and JMod's
parquet has seventeen of the columns its csv has. So the loader intersects what
was asked for with what the file has, and a consumer checks a column arrived
rather than assuming.

**`rename=False` exists for consumers written against the engine's own names.**
This package normalizes every engine onto one schema, which is the right default
and the wrong one for a dashboard whose twenty-five pages say `Retention time`
and `Modified sequence`. Reversing the mapping afterwards is what it was doing.

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

**The cross-repository wiring was removed, deliberately.** For two days this
repository carried a shared index mounted read-only, a self-refreshing working
copy of it, two session hooks, an installer, a membership list, a pasted
`CLAUDE.md` stanza and two workflows called out of a shared repository. It was
built for a network of peers in constant communication. What exists is a
dependency graph with two edges — `streamlit-DO-MS` and `decoder` pin a wheel
from here — and adding an unrelated repository to the account required five steps
in two other repositories, which is the tell that the shape was wrong. So: no
shared repository, nothing mounted, no hooks, and every repository self-contained.
A consumer owns its pin; this repository owns the table of who consumes it,
because that is a fact about this repository rather than a hierarchy. Reading
another repository is `gh repo clone`, which costs a second and cannot go stale.
What survived is what actually fixed the original problem, which was never
communication: tags, wheels, a CHANGELOG, pinned consumers and CI.

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
