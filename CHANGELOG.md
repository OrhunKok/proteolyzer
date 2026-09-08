# Changelog

Notable changes, newest first. Versions are what a consumer pins: the tags this
file names are the ones `pip install` should be pointed at, rather than a commit.

Until 1.0 a minor version may break an interface. What breaks is listed here,
with what to do about it, because two repositories depend on this one -- see the
table in CLAUDE.md -- and the first they knew of the last rename was an
ImportError.

## Unreleased

### Changed

- The history and the measurements have moved out of the source into
  [`docs/notes/`](https://OrhunKok.github.io/proteolyzer/notes/) — what was
  measured off a real export, what a benchmark cost and what its confounds were,
  an assumption that shipped and was wrong. `core.formats`, `core.loader` and
  `core.models` lost ~140 lines of that, so a docstring now says what the thing
  does and points at the note where the *why* is long. `DECISIONS.md` is back to
  a claim per entry, which is what it says it is, and `CLAUDE.md` gained a
  pointer rather than more prose.

  Nothing about behaviour changes. It is here because docstrings are the API
  reference, so the reference reads differently.

### Removed

- `core.io` (`read_frame`, `write_frame`, `frame_exists`) and `core.pipeline`
  (`Stage`, `NullQueue`). They were parquet interchange and stage plumbing for
  the out-of-tree pipelines that no longer exist, and nothing used them:
  not this package, and not either consumer — checked by cloning both and
  grepping, rather than by trusting a search index.

  Dead public surface is not free. It is tested, type-checked, documented, and
  has to be considered for compatibility on every release. This repository's
  own rule is that code lives here when a consumer wants it; the same rule
  removed `LOAD_COLS` in v0.4.0 and refused `quirks.py`.

  **What a consumer has to do.** Nothing, unless you import one of those five
  names — in which case they were 55 and 95 lines, and lifting them into your
  own repository is the whole of the migration. `core.logging` is unaffected
  and stays: it is what every class here logs through.

### Changed

- The README and the docs no longer point at `proteolyzer-cellenone` and
  `proteolyzer-aas`. Neither repository exists — both 404 — and the claim built
  on them was wrong twice over: reading a cellenONE preparation is **in** this
  package (`proteolyzer.cellenone`, since v0.2.0), and the
  amino-acid-substitution pipeline is nowhere. The package description said the
  package "provides the domain pipelines used for ... amino acid substitution
  discovery", which it has not since v0.1.0, and that description ships with the
  wheel.

  Also noted rather than changed: `core.io`, `core.logging` and `core.pipeline`
  were documented as the surface those pipelines used. Nothing in this
  repository or in either consumer uses them today, so the docs now call that
  surface available rather than load-bearing.



### Fixed

- **A float is read exactly, whichever parser reads it.** pandas' own float
  parser stops at about sixteen significant digits, so the q-value
  `0.0007162974636774825` came back `0.0007162974636774` — a relative error of
  1.2e-13, four orders of magnitude past float64's epsilon. It affected the two
  places pandas parses rather than Arrow: the stock CSV parser, which this
  package falls back to when memory is short, and `to_numeric`, which read the
  columns a format writes as text. Both are exact now, so the fallback parser is
  an optimization in reverse rather than a different answer — which is what
  `core.loader` already said it was.

### Changed

- Columns a format writes as text are parsed by Arrow rather than by
  `to_numeric`, which is **~6x quicker** on that step and **1.34x on a whole
  parquet read** of a 170,795-row report. It is a vectorized cast over the
  buffer the values already sit in rather than a walk, and it is the more
  accurate of the two, per above. `to_numeric` remains the fallback for anything
  Arrow will not take.



### Fixed

- **A word is not a gap.** `read_csv` nulls `NA`, `N/A`, `None`, `null` and
  `NULL` by default, and this package took that default — so a text column
  holding one of those words as a *value* came back empty. On a real Spectronaut
  export `EG.InSourceFragmentationClass` is a three-state classification:

  | value | rows | meaning |
  |---|---|---|
  | `Likely Parent` | 1,117 | fragmented in the source |
  | `Likely Child` | 1,146 | a fragment of one, with its parent named |
  | `None` | 168,532 | neither |

  The third state was being erased into "not recorded" for 99% of the report,
  and the same file read as parquet kept it — so the two serializations
  disagreed about data, not just dtype. The readers now null exactly what
  v0.13.0 defined as a gap, which is the machine artefacts (`NaN`, `<NA>`,
  `#N/A`, `1.#IND`, an empty field) and no words at all. One definition, both
  paths.

  **What a consumer has to do.** Check any column you read from a delimited or
  Excel file that uses a *word* for absence. `NA`, `None`, `null`, `NULL` and
  `N/A` now arrive as those strings rather than as gaps — so `.isna()` over such
  a column returns False where it used to return True, and a numeric column that
  spells its holes `NA` (rather than leaving them empty) will come back as text.
  An empty field is still a gap, so a column with ordinary holes in it is
  unaffected; the four search engines whose output is read here all leave theirs
  empty.



### Fixed

- **Text standing for a gap is read as a gap.** A delimited reader turns `NaN`
  and an empty field into missing values itself; parquet stores the string it
  was given. So the same report read the two ways disagreed about which cells
  were empty — and the parquet one lied the worse way round, since `pd.isna`
  says `False` of the string `"NaN"`, so a gap read as data all the way to
  whatever plotted it. On a real Spectronaut export that was two columns,
  one of them empty in every one of its 170,795 rows.

  Only the two spellings that are not words. `NA` and `None` are left as they
  were found: on that export `EG.InSourceFragmentationClass` is `"None"` in
  168,532 rows of 170,795 with a real class in the rest, which makes it a
  category and not an absence. A text reader nulls both by default, so the two
  serializations still disagree about those two — agreeing means deciding
  pandas' default is right, which is a wider change than this one.

- Detection no longer raises where a file's columns match two formats. Two
  blocks naming one file is this package's config contradicting itself and
  still raises; two signatures matching is a fact about somebody's file — a
  report joined to another engine's table, saved, read back — and that now
  warns and reads as `Unknown`, which is what it did before v0.11.0.

### Changed

- The guide has a page on [DIA isolation windows](guide/isolation.md).
  `envelope_split` and `envelope_room` have been public since v0.3.0 and v0.7.0
  and appeared only in the generated API reference, so a reader of the docs had
  no way to learn the question could be asked. It says what each answers, what
  the caller has to supply and why, and that `Window` is a position.
- `CHANGELOG.md` said three repositories depend on this one; there are two.
  The quickstart in `README.md` and the docs still called recognition a matter
  of the file name, which it has not been since v0.11.0.



### Added

- **A table too narrow to sign by an overlap of names is claimed by its whole
  shape.** v0.11.0 recognized four of DIA-NN's five tables by their columns and
  said the `xic` export could not be: `pr`, `feature`, `rt`, `value` are four
  words that belong to nobody, and one of them is JMod's. That was giving up too
  early — what separates them is not *which* of those names is present but that
  **all** of them are and there is almost nothing else, four columns against
  JMod's thirty-four.

  A format block may now carry `NARROW_SIGNATURES`, and both halves have to
  hold: every declared column present, and the table about as narrow as
  declared. Either half separates this pair on its own, which is the argument
  for both — requiring the whole of a small schema is what works against another
  engine sharing a name, and the width is the second lock against a frame that
  carries those four words among two hundred others. Shape is asked before the
  overlap, being the more specific claim, and both are still behind the name.

  With it, the invariant test covers **every** real table of every engine with
  no carve-out: twelve tables, one claimant each.

  **What a consumer has to do.** Nothing; this reaches files that came back
  `Unknown`. It removes the one limitation v0.11.0 shipped with.

## v0.11.0

### Added

- **Every format is recognized by its columns, not just by its file name.**
  A file name is a convention people depart from — you rename what you
  downloaded, and a report called `2026-08-27_experiment_three.parquet` used to
  come back `Unknown`, unrenamed, with a warning. DIA-NN, MaxQuant, JMod and
  FragPipe now each carry a `COLUMN_SIGNATURE` alongside Spectronaut's, matched
  against the file's own columns when no format claims the name.

  The name is still asked first, so a file called what its engine calls it is
  claimed without being opened — cheaper, and identical to every release before
  this one.

  Bounded by an invariant test that walks **every real table of every engine**
  and asserts exactly one format claims it. That matters because detection
  refuses a file two formats claim: a signature reaching another engine's table
  would not mis-read one file, it would make both unreadable. Four names are
  shared between engines and appear in no signature — `Charge` and `Intensity`
  (MaxQuant and FragPipe), `PEP` (MaxQuant and DIA-NN), `rt` (JMod and DIA-NN's
  XIC export).

  **What a consumer has to do.** Nothing. This only reaches files that came back
  `Unknown`; anything recognized before is recognized the same way. Two limits
  worth knowing: a name that *does* match still wins, so a table deliberately
  renamed to another engine's filename is read as that engine; and DIA-NN's
  `xic` export cannot be signed — its four columns are `pr`, `feature`, `rt`,
  `value`, one of which is JMod's — so it stays recognized by name alone.

## v0.10.0

### Fixed

- **A column a search engine writes as text comes back with its dtype.**
  Spectronaut serializes some numbers and flags as strings, inconsistently
  within the one file: `EG.Qvalue` arrives as `'1.99e-13'` while `PG.Qvalue`
  beside it is a double, and `PEP.IsProteotypic` is `'False'` while
  `EG.IsDecoy` is a real boolean. So the canonical `Q.Value` was a string, and

  ```python
  frame["Q.Value"] < 0.01
  TypeError: Invalid comparison between dtype=str and float
  ```

  which is the first thing anyone does with a report. `Q.Value`,
  `Missed.Cleavages` and the channel q-values now come back numeric and
  `Proteotypic` boolean, matching what the same report read from text gives.

  The columns are a list on the format block rather than a rule, because the
  rule gets it wrong: `FG.XICDBID` is a database key whose every value parses
  as a number, and turning an identifier into an integer is quiet damage.
  Conversion is **all or nothing** per column -- a column where any value is
  not a number is left exactly as it arrived and logged, rather than being
  handed back numeric and shorter by however many values nobody was told
  about. A column an export already stored properly is not touched, and the
  text a file uses for a gap (`'NaN'`, which `pd.isna` calls False) becomes a
  real gap.

  Applied before the rename, on the file's own column names, so it reaches a
  caller reading with `rename=False` -- `streamlit-DO-MS` reads every format
  that way, and a q-value that cannot be compared to a float is no more use
  under one name than another.

## v0.9.0

### Fixed

- **The Spectronaut rename mapping did not fire on a real parquet export.**
  v0.8.0 added `.parquet` to the format and kept the tab-separated export's
  column names, on the assumption that the two spell them the same way. They do
  not: the parquet export writes `R_FileName` where the text one writes
  `R.FileName`, and turns the space in `PG.Cscore (Run-Wise)` into an underscore
  as well — a dot being a path separator in a nested parquet schema. So every
  one of the seventeen mapped columns missed, and a report came back under the
  file's own names with nothing renamed and no `Precursor.Id` built.

  Both spellings are mapped now, derived from one list so they cannot drift.
  Checked against a real 170,795-row export: all seventeen match.

### Added

- **A Spectronaut report is recognized by its columns, not its file name.**
  Spectronaut has no default output name — whoever runs the analysis names the
  export, so `GluC-30min.parquet` is as real a report as any, and the
  `..._Report` pattern v0.6.0 leaned on was a convention rather than a rule.
  A format block may now carry a `COLUMN_SIGNATURE`, and where no block claims
  the name, the file's own columns are read and matched against it: a parquet
  footer or one header line, so the usual case pays nothing.

  This is the call the cellenONE reader already makes for the same reason —
  which file is which is worked out from the file, because names are unreliable.
  Only Spectronaut carries a signature, and two of its columns have to match, so
  looking inside cannot start claiming another engine's output.

  **What a consumer has to do.** Nothing. A file recognized before is recognized
  the same way, by name, without being opened; this only reaches files that came
  back `Unknown`.

## v0.8.0

### Added

- **Spectronaut reports are read as parquet as well as tab-separated text**,
  parquet being what Spectronaut writes by default. Same names, same rename
  mapping, same built `Precursor.Id`; the extension decides which reader runs
  and a caller does nothing to pick. There is a test asserting the same report
  written both ways comes back as the same frame, values and dtypes.

  v0.6.0 claimed `.tsv` alone because the export the format was measured from
  was tab separated, which made the *default* export unrecognized: it read as
  `Unknown`, so the frame came back under Spectronaut's own column names with a
  warning rather than onto the canonical schema.

  **What a consumer has to do.** Nothing, unless it was reading a
  `..._Report.parquet` as an unrecognized file, which now comes back renamed —
  `rename=False` keeps the file's own names, as `streamlit-DO-MS` does.

  Note that DIA-NN also claims `.parquet`, so its `report.parquet` and a bare
  Spectronaut `Report.parquet` are now one capital letter apart on the same
  extension. Detection is case-sensitive and refuses a file two formats claim,
  so this is pinned by a test rather than left to chance.

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
