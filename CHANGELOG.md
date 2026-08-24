# Changelog

Notable changes, newest first. Versions are what a consumer pins: the tags this
file names are the ones `pip install` should be pointed at, rather than a commit.

Until 1.0 a minor version may break an interface. What breaks is listed here,
with what to do about it, because three repositories depend on this one and the
first they knew of the last rename was an ImportError.

## Unreleased

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
