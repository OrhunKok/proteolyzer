# Changelog

Notable changes, newest first. Versions are what a consumer pins: the tags this
file names are the ones `pip install` should be pointed at, rather than a commit.

Until 1.0 a minor version may break an interface. What breaks is listed here,
with what to do about it, because three repositories depend on this one and the
first they knew of the last rename was an ImportError.

## v0.7.0

### Added

- `core.envelope_room(report, windows)`: which window isolated each precursor
  and how much room it left the envelope, rather than only the verdict
  `envelope_split` already gave. Returns a frame of `Window` and `Room` (m/z
  from M+2 to that window's upper edge, negative where the envelope was split,
  NaN where no window isolated it). `envelope_split` is now a thin wrapper over
  it — the sign of `Room` is exactly its verdict — so `streamlit-DO-MS` no
  longer has to re-derive the window index from `ISOTOPE_STEP` and
  `ENVELOPE_ISOTOPES` to say which window a plot should widen.

  **`Window` is a position, so it is `windows.iloc[w]`, not `windows.loc[w]`.**
  The two are the same only while the design still has its original index, and
  a caller that filtered its windows has one that does not — where `.loc` is
  then a `KeyError`, or the wrong window with nothing said. It is -1 where no
  window isolated the precursor, which is a sentinel a position can carry and a
  label could not.

### Fixed

- `envelope_split` (and now `envelope_room`) deduped precursors on ion
  mobility whenever the report carried an `IM` column, even for an m/z-only
  window design that never reads it. Mobility is measured per identification,
  so as a dedupe key it barely dedupes at all — on a measured 1.04M-row
  report, 1,040,711 unique rows against 169,697 without it, a 9x slowdown for
  an identical answer. The key now only includes `IM` when the design itself
  carries `WINDOW_MOBILITY_COLUMNS`.

  What you get back depends on how many identifications your report holds per
  precursor, since that is the whole of what the dedupe was failing to do: 9x on
  that report, and 2.3x measured independently on 300k rows over 5k precursors.
  It bites any m/z-only design read alongside a report that has an `IM` column,
  which is every `--export-windows` file off a timsTOF.

## v0.6.0

### Added

- **Spectronaut is the fifth format read.** `Config().Spectronaut` describes the
  long-format report — one row a precursor a run — and `pz.read` recognizes,
  reads and renames one like any other engine's output.

  ```python
  pz.read("20260901_164751_GluC_30min_Report.tsv")            # onto our names
  pz.read(upload, INPUT_TYPE="Spectronaut", cols_to_load={...}, rename=False)
  ```

  Three things about this format that are not true of the other four:

  **The name is stamped, so it is matched by pattern.** Spectronaut writes
  `<date>_<time>_<analysis>_Report.tsv` and no fixed name could match that, so a
  block may now carry `FILE_PATTERNS` beside `FILES` — a regex over the stem,
  matched in full. Only Spectronaut has one. It is matched case-sensitively,
  because DIA-NN's `report.tsv` is one capital letter from a bare `Report.tsv`
  and a file two blocks claim is refused rather than guessed at.

  **The report carries no precursor identifier**, `EG.PrecursorId` not being in
  every export, so `Precursor.Id` is built out of `EG.ModifiedSequence` and
  `FG.Charge` as the file is read — which is where it has to happen, since
  `process()` asks for that column before any of its own steps run. A block says
  what to build and out of what in `BUILT_COLS`, and `Data.built_cols` is empty
  under `rename=False`: those names are the core's vocabulary, and a caller
  keeping the file's own has not asked for it.

  **A quantity is not on the scale another engine's is** — 2.54 and 400,000 in
  the one column — so `round_large_floats` would take a fifth off the low end of
  it. It is off by default, as it is for every format, and there is now a test
  that says so about this one.

  Written from a measured export: 13 runs, 173,443 rows, 174 MB, tab separated,
  CRLF, 78 columns. What one holds is configurable column by column, so the
  intersection `cols_to_load` already takes is load-bearing here rather than
  convenient — a list written against one lab's export names columns another's
  does not have.

  **What a consumer has to do.** Nothing, unless it was reading a
  `..._Report.tsv` as an unrecognized file: that now comes back under the
  canonical names rather than the file's own. Pass `rename=False` to keep them.
  `streamlit-DO-MS` already reads every format that way.

- `cellenone.CoordinatesMapping.map_data()` carries `ImageFile`, and one column
  per further imaging channel — `ImageFile.Green` beside it. cellenONE
  photographs each cell it prints and names the file in the geoprops table; that
  name was being dropped, and it is the only link from a cell to its picture.
  Unwrapped from the `=HYPERLINK("...")` a spreadsheet writes, so a consumer gets
  a file name relative to the run folder rather than a formula. Missing for a cell
  that was never printed, which is most of them.

## v0.5.0

### Added

- `cellenone.CoordinatesMapping.map_data()` carries a `Pickup.Source` column: the
  pickup file (or Run ID) whose dispensing placed each cell. `Plate` is the
  destination plate's position in the run rather than its identity — one plate
  is mounted at a time, so it reads 1 for every real pickup — so a preparation
  spanning several plates had no way to tell `A1` on one from `A1` on another.
  `Pickup.Source` does.

### Changed

- `openpyxl` is a dependency rather than absent. `.xlsx` is a routed extension in
  `core.loader`, so a caller handed a spreadsheet was told the format is supported
  and then given an `ImportError` from pandas. It is ~250 KB of pure Python and
  pandas imports it inside `read_excel`, so importing proteolyzer still does not
  pay for it — `tests/test_package_boundaries.py` keeps that true.

### Fixed

- The droplet-table fixture writes short rows, as the instrument does. Nothing
  changes for a caller; what changes is that the suite would now catch the defect
  released in v0.2.2 — a reordered geoprops export read as a second table of
  cells — instead of needing a real run directory to find it.

## v0.4.0

### Changed — breaking

- **The core no longer decides which columns you get.** `LOAD_COLS` is gone from
  every format block, and a file is read whole unless the caller names columns.
  `Data.cols_to_load` is where a project states its own list.

  Which columns matter is a fact about the project reading the file, not about
  the file: a dashboard plots the m/z and injection time a pipeline never looks
  at, and a pipeline wants quantities no panel shows. The 260 names across 29
  tables that lived here were contributed by one consumer, which had already
  overridden two of them because "the core's own subset is a pipeline's rather
  than a dashboard's". A shared list that every sharer overrides is not shared.

  **What a consumer has to do.** If you passed `cols_to_load`, nothing changes —
  you already stated your own list. If you passed `extra_cols_to_load` and
  relied on it meaning *the core's subset plus mine*, you will now get every
  column: correct, and wider. Move that list to `cols_to_load` to get a narrow
  frame back. To recover exactly what you had, install v0.3.0 and dump it:

  ```python
  from proteolyzer.core.formats import Config
  {t: sorted(c) if c else None
   for t, c in Config().DIANN.LOAD_COLS.items()}   # or MaxQuant, JMod, FragPipe
  ```

- `extra_cols_to_load` on its own no longer narrows anything, there being no base
  subset for it to be extra to. It widens `cols_to_load` where both are given.

### Removed

- `Config().<engine>.LOAD_COLS`. Format recognition (`FILES`, `FILE_EXTENSIONS`),
  the canonical rename mapping and the categorical exclusions all stay — those
  are facts about the format, true for everyone who reads it.

## v0.3.0

### Added

- `core.envelope_split` says whether a precursor's isotopic envelope was isolated
  whole by a DIA window design, or split across the edge of one. Where part of the
  envelope is fragmented in another window, the MS2 spectrum covers less of the
  signal than the MS1 quantitation was measured over — a property of the method
  and the precursor rather than of the sample, and worth ruling out before a
  quantitative difference is read as biology.

  Moved from streamlit-DO-MS, which had the only copy and could ask the question
  only from inside a dashboard. Verified against the original over 400 random
  window designs before the copy there is dropped.

  Note that `Precursor.Mz` is not in the core's DIA-NN column subset, that subset
  being a pipeline's rather than a dashboard's, so the caller supplies the frame:
  one without the columns is answered with nothing rather than a guess.

- `core.jaccard_index`, how far two masks agree — what they share over what either
  has. Also from streamlit-DO-MS, where it compares the precursors two labelling
  channels or two runs identified. `nan` for two empty masks, where the original
  divided by zero.

## v0.2.2

### Fixed

- A cellenONE folder yields one table of cells, not two. The instrument writes
  two tables with the same fifty-two columns — the cells it printed, and every
  droplet it detected on the way — and they are told apart by what fraction of
  their rows say which target, field and position a cell went to. That fraction
  was computed after short rows had been discarded, and a droplet that was never
  placed is exactly a short row: of the 6797 rows one real geoprops table had,
  6770 were dropped before the ratio was taken, which then came to 1 for both
  tables. So a reordered export was read as a second table of cells in every
  folder that had one. Only the half line a capped read ends on is dropped now.
  A consumer reading a run with reordered exports gets fewer cells than before,
  and every cell it loses is one the instrument never printed.

## v0.2.1

### Fixed

- A cellenONE run directory can be read twice. Every parser now starts the file
  at its start rather than wherever the last reader left it: an upload is read to
  work out which step each file belongs to and read again to parse it, and a
  consumed buffer parses as an empty file — which arrives as 'No columns to parse
  from file'. It is the second read that matters, because correcting a step and
  reading the run again is what the correction is for.
- `map_data()` parses the droplets if they have not been parsed already, rather
  than raising `AttributeError` about an attribute a caller has no reason to know
  exists.

## v0.2.0

### Added

- `proteolyzer.cellenone` reads a cellenONE run directory: which file belongs to
  which step of the preparation, one row per cell with where it landed and what
  was dispensed onto it, and the chamber readings stamped with the step they were
  taken during. See the [guide](guide/cellenone.md).

  It came from streamlit-DO-MS, which had the more developed of two copies of it:
  the same class name and method names existed in a second repository, both had
  independently fixed the same imaging-channel bug, and neither knew about the
  other. The two are now one, and it lives where the projects that read this
  instrument's output can share it.

- A fluorescence channel is kept rather than discarded. cellenONE writes a
  geoprops row per cell *per imaging channel*, so counting rows doubles the cell
  count and halves every geometry average. Both copies had noticed; one filtered
  to the transmission row, which fixes the count and throws the fluorescence
  away. Here the transmission row is the cell and every other channel joins onto
  it under its own name — `Diameter.Green` beside `Diameter` — because in a
  sorting experiment the fluorescence is what the experiment was for. A
  preparation imaged in transmission alone comes out exactly as before.

- The cells table is recognised by its columns whatever channel it was imaged
  in, so a fluorescence-only export is no longer invisible.

### Note for consumers

The subpackage is imported on first use, like `plots` and `unimod`, so a core
install does not pay for it. It needs nothing that the core does not.

## v0.1.0

The first tagged release. It exists because the packages that depend on this one
were pinning commits, or nothing at all: `decoder` broke at import time on a
rename with no version boundary to notice (decoder#3), and `streamlit-DO-MS` was
pinned to a 40-character SHA.

The state of the interface as of this tag:

### The core

- `proteolyzer.core` — file recognition, reading, dtype narrowing, derived
  columns, labelling, and pivoting to a quantitative matrix. Was
  `proteolyzer.utils`.
- Domain constants are `proteolyzer.reference`; they were `proteolyzer.config`,
  which no longer exists. `pz.config.Codons.Standard` is
  `proteolyzer.reference.CODON_TABLE`, a `Mapping[str, str]`.
- `Data.load()` returns a `Report`, not a DataFrame. `Report.frame` is the
  DataFrame; the wrapper carries where the data came from and what was done to
  it. `Report` is deliberately not a DataFrame subclass.
- Four search engines are recognised: DIA-NN, MaxQuant, JMod and FragPipe. A
  format is a block on `core.formats.Config`, and detection walks whatever
  blocks are there.
- `Data` takes `cols_to_load` to replace the configured column subset,
  `extra_cols_to_load` to add to it, and `rename=False` to keep the file's own
  column names.
- `narrow(frame)` and `Narrower` narrow a frame's dtypes without the pipeline
  around them, for a caller that derives its own columns.
- `Report.matrix(...)` pivots to a quantitative matrix; missingness and group
  normalisation live there.

### What moved out

- The cellenONE and amino-acid-substitution pipelines left this repository for
  their own. Note for anyone tracking master: cellenONE is on its way back in,
  because two projects were maintaining separate copies of it.

### Fixed, for anyone who pinned a commit before this tag

- `msScans` was listed in MaxQuant's load columns but not in its files, so a
  `msScans.txt` was not recognised as MaxQuant output and came back as one
  undelimited column.
- `extra_cols_to_load` was ignored for a file with no configured subset, so
  asking for two columns of it read all of them.
- `allPeptides`, `msScans` and `msmsScans` had no column subsets, and a file
  with no subset is read whole.
- The run column in every MaxQuant table except `evidence` is `Raw file`, which
  the rename mapping did not reach.
- A source archive can be built: `setuptools_scm` has a fallback version, so
  `pip install <archive>.tar.gz` no longer fails for want of a `.git`.

### Depending on this

Pin the wheel attached to the release, not the source archive:

```
proteolyzer @ https://github.com/OrhunKok/proteolyzer/releases/download/v0.1.0/proteolyzer-0.1.0-py3-none-any.whl
```

A GitHub source archive carries no `.git`, so setuptools-scm falls back and the
installed version reads `0.0.0` whichever tag it came from — pip cannot then
tell one release from another. The wheel carries the version in its metadata.
`.github/workflows/release.yml` builds and attaches both on every tag.

`make test-downstream` runs the suites of the repositories that depend on this
one against the working tree, when they are checked out under `downstream/`.
Run it before tagging.
