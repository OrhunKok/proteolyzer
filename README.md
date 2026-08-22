![Proteolyzer](./proteolyzer_logo_with_name.svg)

---------

[![CI](https://github.com/OrhunKok/proteolyzer/actions/workflows/ci.yml/badge.svg)](https://github.com/OrhunKok/proteolyzer/actions/workflows/ci.yml)

# Proteolyzer

Proteolyzer is a Python package for processing, analyzing, and visualizing
proteomics data. It reads DIA-NN and MaxQuant output, normalizes it into a
consistent shape, and provides the domain pipelines used for single-cell
sample preparation and amino acid substitution discovery.

**Documentation:** <https://OrhunKok.github.io/proteolyzer/>

## Features

- **Data loading** — Parquet, TSV/CSV, Excel and plaintext, with per-format
  column subsetting so only the columns you need are read, through pyarrow's
  multithreaded parser where it can be used.
- **Data processing** — dtype narrowing, derived columns, missed-cleavage
  flags, per-run identification counts, and automatic detection of labelled
  (mTRAQ/SILAC/TMT) precursors.
- **Matrix transformation** — pivot to a quantitative matrix, measure
  missingness, and normalize within column groups.
- **Plotting** — publication-styled relational plots (volcano plots with
  significance and effect-size thresholds, point labelling).
- **UniMod plugin** — query the full UniMod database with SQL; the database is
  built on first use rather than shipped.

## Installation

Python 3.14 or newer.

```bash
git clone https://github.com/OrhunKok/proteolyzer
cd proteolyzer
pip install .                 # core
pip install '.[unimod]'       # only to *build* the UniMod cache; querying needs nothing
pip install -e '.[dev]'       # editable install with test and lint tooling
```

To depend on it, pin a tag rather than a commit or a branch — see
[CHANGELOG.md](./CHANGELOG.md), which says what each one broke:

```
proteolyzer @ https://github.com/OrhunKok/proteolyzer/archive/refs/tags/v0.1.0.tar.gz
```

Notes on the optional extras:

- The plotting module defaults to the `science`
  ([scienceplots](https://github.com/garrettj403/SciencePlots)) theme, which
  renders text with LaTeX. Without a LaTeX installation, pass another theme:
  `VolcanoPlot(..., theme="default")`.

## Quickstart

### Core: load, process, pivot

```python
import numpy as np
import proteolyzer as pz

report = pz.read("report.parquet")   # DIA-NN report, recognized by name
processed = report.process()         # dtypes, derived columns, labelling info

processed.runs                       # {'run1', 'run2', ...}
processed.n_identifications          # distinct precursors
processed.summary()                  # precursors/peptides/proteins per run
processed.processing.label_free      # what processing found
processed.frame                      # the pandas DataFrame

matrix = (
    processed.matrix("Ms1.Area", index=["Precursor.Id"], columns=["Run"])
    .normalize_matrix(within_groups=["Run"], agg_func=np.nansum)
    .matrix
)
```

`read` recognizes known DIA-NN/MaxQuant file names and reads only the columns
configured for them. Use `load_all_columns=True` for everything, or
`extra_cols_to_load=[...]` to add to the default set.

A `Report` is a frozen wrapper around three things: the `frame`, the `source`
it was read from, and the `processing` that produced it (`None` until
`process()` runs). It deliberately is *not* a DataFrame subclass — pandas
returns plain frames from most operations, so a subclass loses its metadata
the first time anyone slices it. `len()`, `report["Run"]` and `.columns` work
for interactive use; everything else goes through `.frame`, which is always a
plain `DataFrame`.

#### Memory and dtypes

`process()` works to make a report cheaper to hold, which on a real one cuts
it by a third to a half:

- Float columns that already hold whole numbers become the smallest exact
  integer dtype. Columns with a real fractional part keep it; pass
  `process(round_large_floats=True)` to also round high-magnitude columns
  (intensities) to integers, discarding their fractional part.
- Remaining `float64` columns narrow to `float32`, which is what DIA-NN
  stores in its own parquet output — so a report read from the TSV ends up
  costing the same as the same report read from parquet, instead of 44% more.
  The worst relative error measured across every float column of a real
  report is 6e-8, float32's epsilon. A column is left alone if any value
  falls outside float32's range, where narrowing would turn it into an
  infinity or a zero. Pass `process(narrow_floats=False)` to keep double
  precision — worth doing if you subtract nearly equal values from the frame
  yourself, since that amplifies the rounding.
- Text columns become categorical where measuring both says that saves at
  least 20% of the column. Protein and gene columns typically halve;
  identifier columns, which barely repeat, are left as they are.

`report.memory()` gives the per-column breakdown, largest first, if a report
is bigger than you expect.

Reading is chosen the same way. The pyarrow parser is several times faster
than the stock one but holds its Arrow table and the frame at once, so it is
used only when the file is small enough for that against the memory
available; past that the stock parser reads it with a bounded footprint.

### Logging

Progress is reported on the `proteolyzer` logger, which is configured on import
and never touches the root logger:

```python
import logging
import proteolyzer as pz

pz.configure_logging(level=logging.DEBUG)            # more detail
logging.getLogger("proteolyzer").setLevel("WARNING")  # quieter
```

### Reference data

Masses, the genetic code and protease rules live in one place, as immutable
mappings read on first use:

```python
from proteolyzer import reference

reference.amino_acid_masses()["A"]        # 71.03711...
reference.protease("Trypsin").allowed_counts
reference.CODON_TABLE["AUG"]              # "M"
```

They are exported from UniMod. For anything the two exported tables cannot
answer, query the database itself:

```python
from proteolyzer import unimod

unimod.tables()                                    # built on first use
unimod.table("modifications")
unimod.query(
    """
    SELECT m.full_name, s.one_letter, m.mono_mass
    FROM specificity AS s
    JOIN modifications AS m ON s.mod_key = m.record_id
    WHERE m.mono_mass BETWEEN ? AND ?
    """,
    (15.9, 16.1),
)
```

The database is 5 MB against 221 KB of CSVs, so it is not shipped: it is built
once into `$PROTEOLYZER_CACHE_DIR` (else `$XDG_CACHE_HOME`, else `~/.cache`)
and reused. Querying uses only the standard library; building downloads the
UniMod XML and needs `pip install '.[unimod]'`.

Regenerating the bundled CSVs after a new UniMod release is a maintainer step:

```bash
python -m proteolyzer.unimod build
python -m proteolyzer.unimod export \
    --mods-output src/proteolyzer/resources/unimod_modifications.csv \
    --aa-output   src/proteolyzer/resources/unimod_amino_acids.csv
```

## Project layout

```
src/proteolyzer/
    reference.py      domain constants: masses, codons, proteases
    core/             the base suite, imported eagerly
        formats.py    the input formats recognized, and their columns
        loader.py     read a described source into memory
        models.py     Data (a described input) and Report (what comes back)
        processor.py  dtypes, derived columns, labelling information
        matrix.py     pivot to a quantitative matrix
        io.py         parquet interchange between pipeline stages
        pipeline.py   stage plumbing: parameters, progress, provenance
    plots/            optional: plotting base class and relational plots
    unimod/           optional: UniMod SQL queries, built on demand
tests/                pytest suite
examples/             runnable notebooks per module
docs/                 documentation pages; the API reference is generated at build
```

The core and `reference` are imported with the package. `plots` and `unimod`
load on first access, so a missing extra only fails for the module that needs
it — `tests/test_package_boundaries.py` enforces that, including that importing
proteolyzer does not pull in matplotlib.

## Domain pipelines

Pipelines for a particular instrument or assay live in their own repositories,
built on this core, so that neither their dependencies nor their release
cadence lands on everyone installing it:

| repository | what it does |
|---|---|
| [proteolyzer-cellenone](https://github.com/OrhunKok/proteolyzer-cellenone) | maps single cells prepared on a cellenONE to well positions, and flags well/label clashes |
| [proteolyzer-aas](https://github.com/OrhunKok/proteolyzer-aas) | discovery of amino acid substitutions ([paper](https://decode.slavovlab.net/)) |

Each uses `core.io` for the parquet interchange, `core.logging` for the logger,
and `core.pipeline` for stage plumbing and provenance. That surface is what to
keep stable.

## Development

```bash
pip install -e '.[dev]'
pre-commit install                  # run the lint hooks on commit
pytest                              # test suite
ruff check .                        # lint
ruff format src tests               # format
mypy                                # type check (scope in pyproject.toml)
mkdocs serve                        # documentation, with live reload
```

`make help` lists the same targets.

Tests that need an optional dependency skip themselves, so the suite runs on a
core install. Run it from the repository root.

## Contributing

Contributions are welcome! Please fork the repository, make your changes, and
submit a pull request. Ensure `ruff check .` and `pytest` pass, and add tests
for new behaviour.
