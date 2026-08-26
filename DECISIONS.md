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

**`Plate` is a position, not an identity, so `map_data()` also carries
`Pickup.Source`.** One destination plate is mounted at a time, so the column
the instrument writes reads 1 for every real pickup — a prep spanning several
plates writes 1 for each of them, and their wells then collide: `A1` from
plate one and `A1` from plate two become one key. The file that did the
dispensing is what tells them apart, and `log_parse` already stamps every row
with it to settle a redo against the attempt it replaced; `Pickup.Source` is
that same value, kept instead of dropped, carried through the nearest-cell
match onto the geoprops row it placed.

**The image of a printed cell is carried, and unwrapped on the way.** cellenONE
photographs each cell it prints and names the file in the geoprops table, which is
the only link from a row to a picture -- nothing else in an export carries one, and
it used to be dropped on the way through `xls_parse`. It arrives as a spreadsheet
formula rather than a path, `=HYPERLINK("26_Printed_T1_F1_R9_C1_(A-1)_Trans...")`,
so it is unwrapped to the bare name, which is relative to the run folder.

Only a *printed* cell is photographed: on the run this was written from, 1,260 of
160,091 rows carried one, being 630 cells in two channels. So a caller gets the
transmission view and the viability view of a cell that was kept, and nothing at
all for the ones that were not -- which is worth knowing before anyone plans to
learn viability from these images, because the dye is read in exclude mode and the
excluded cells are never photographed.

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

**The agent opens pull requests as a GitHub App, not as a person.** A pull
request opened with `GITHUB_TOKEN` gets no checks: Actions does not fire
workflows for events that token creates. `ci.yml` triggers on `pull_request` and
still never ran on one the agent opened, so the gate ran inside the agent's own
job and the pull request showed a blank check list — work that had been tested,
with no evidence of it where a reviewer looks.

A personal access token fixes that and costs two things worth more. It is scoped
to everything the account owns rather than to one repository, and it makes the
agent a `User`. The guard on the agent job is `github.event.issue.user.type !=
'Bot'`, and the model here is that a change belonging elsewhere is filed as an
issue in that repository — so with a PAT, an issue this agent files is picked up
by that repository's agent as though a person had filed it. That is the loop the
guard exists to stop, and a PAT opens it while looking like a fix.

An App installation token is minted per run, dies with it, is scoped to the
repository, and is still a `Bot`. With no App configured the step fails and the
run continues on `GITHUB_TOKEN`: the issue is still worked and the pull request
still opened, and the only thing lost is the checks on it.

**The App is installed per repository, and the fallback says so.** "All
repositories" would mean a project added later works without anyone remembering,
which is the whole of the argument for it. Against: the App's reach is what a
leaked private key reaches, and the account holds repositories that will never
have an agent, including forks of other people's work. Installing on the three
that do keeps the key's blast radius to the three it is already stored in.

That leaves one hazard, and it is the one this step exists to prevent arriving by
another door: a fourth project, the App not installed, the token step failing,
the run carrying on, and a pull request with no checks that looks exactly like a
pull request that passed. So the fallback annotates the run. Silence was the
defect; falling back is not.

## What may span repositories

**A credential may be shared. Content may not.** What was torn out on 2026-08-23
was a repository the others ran code from — an installer, a capability index, two
session hooks, a reusable workflow. Its defect was never that it was shared. It
was that it was *content which had to be kept current*, and misled silently when
it was not: the index once named a member living on a branch that had been
deleted.

One thing added since spans all three repositories: the `personal-claude-agent`
GitHub App, so their pull requests carry checks. That is a shared credential, and
it is worth naming as the exception rather than letting it pass unremarked. It is
allowed because it has no content — nothing to pull, nothing to go stale, no
version to be behind. Delete it and every repository still works: the token step
fails, the run carries on with `GITHUB_TOKEN`, and an annotation says so. The test
for the next thing that wants to span repositories is that question. *Can it be
out of date?* A credential cannot. A column list, a capability index and a
reusable workflow all can, and every one of them was.

**Code moves here only when a second consumer wants it.** `envelope_split` came
because it reads the canonical schema and any DIA pipeline wants it. `quirks.py`
was proposed the same day and refused: its helpers read each engine's *raw*
column names, which exist only because streamlit-DO-MS loads with `rename=False`,
and decoder mentions no FragPipe, MaxQuant or JMod anywhere. One consumer means
"written once, upstream" buys a release to coordinate and nothing else — the same
argument that removed `LOAD_COLS` in v0.4.0. That call is closer than it looks
and is worth making out loud each time, because the pull toward putting things
here is exactly what built the last network.

**An issue filed elsewhere is not the bus.** The old one was two labels and a
workflow calling a shared worker, and it needed all three repositories to agree
on it. Opening an issue in the repository a change belongs to is GitHub, not
infrastructure: nothing here reaches out, nothing there depends on this. Pushing
a branch to another repository *is* reaching out, and was done once, on
2026-08-24, to give `streamlit-DO-MS` and `decoder` the App token step — because
their own agents are refused write access to `.github/workflows/` and could not
do it themselves. That is the shape of the exception: a person, once, where the
sanctioned path is closed by design.
