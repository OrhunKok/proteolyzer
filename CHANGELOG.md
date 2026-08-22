# Changelog

Notable changes, newest first. Versions are what a consumer pins: the tags this
file names are the ones `pip install` should be pointed at, rather than a commit.

Until 1.0 a minor version may break an interface. What breaks is listed here,
with what to do about it, because three repositories depend on this one and the
first they knew of the last rename was an ImportError.

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

```
proteolyzer @ https://github.com/OrhunKok/proteolyzer/archive/refs/tags/v0.1.0.tar.gz
```

`make test-downstream` runs the suites of the repositories that depend on this
one against the working tree, when they are checked out under `downstream/`.
Run it before tagging.
