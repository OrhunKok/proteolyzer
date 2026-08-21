![Proteolyzer](./proteolyzer_logo_with_name.svg)

---------

[![CI](https://github.com/OrhunKok/proteolyzer/actions/workflows/ci.yml/badge.svg)](https://github.com/OrhunKok/proteolyzer/actions/workflows/ci.yml)

# Proteolyzer

Proteolyzer is a Python package for processing, analyzing, and visualizing
proteomics data. It reads DIA-NN and MaxQuant output, normalizes it into a
consistent shape, and provides the domain pipelines used for single-cell
sample preparation and amino acid substitution discovery.

## Features

- **Data loading** — Parquet, TSV/CSV, Excel and plaintext, with per-format
  column subsetting so only the columns you need are read.
- **Data processing** — dtype narrowing, derived columns, missed-cleavage
  flags, and automatic detection of labelled (mTRAQ/SILAC/TMT) precursors.
- **Matrix transformation** — pivot to a quantitative matrix, report
  missingness, and normalize within column groups.
- **Plotting** — publication-styled relational plots (volcano plots with
  significance and effect-size thresholds, point labelling).
- **cellenONE module** — maps single cells prepared on a cellenONE to well
  positions, and flags well/label clashes.
- **UniMod plugin** — query the full UniMod database with SQL; the database is
  built on first use rather than shipped.
- **Alternate RNA decoding (AAS) module** — the pipeline used for discovery of
  amino acid substitutions and PTMs ([paper](https://decode.slavovlab.net/)).

## Installation

Python 3.14 or newer.

```bash
git clone https://github.com/OrhunKok/proteolyzer
cd proteolyzer
pip install .                 # core
pip install '.[aas]'          # amino acid substitution pipeline
pip install '.[unimod]'       # only to *build* the UniMod cache; querying needs nothing
pip install -e '.[dev]'       # editable install with test and lint tooling
```

Notes on the optional extras:

- `aas` pulls `quickdna` from a pinned wheel because PyPI has no build for the
  supported interpreter. It is only needed by
  `proteolyzer.aas.translation.FrameTranslator`; the rest of the AAS pipeline
  imports and runs without it.
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

`process()` narrows float columns that already hold whole numbers to the
smallest exact integer dtype. Columns with a real fractional part are left
alone; pass `process(ROUND_LARGE_FLOATS=True)` to also round high-magnitude
columns (intensities) to integers, which discards their fractional part.

### Logging

Progress is reported on the `proteolyzer` logger, which is configured on import
and never touches the root logger:

```python
import logging
import proteolyzer as pz

pz.configure_logging(level=logging.DEBUG)            # more detail
logging.getLogger("proteolyzer").setLevel("WARNING")  # quieter
```

### cellenONE: map cells to wells

```python
import proteolyzer as pz

mapper = pz.cellenone.CoordinatesMapping(
    root_dir="cellenone_files", label_type="mTRAQ", plex=2
)
metadata = mapper.map_data()   # one row per isolated cell, with pickup well and label
stats = mapper.map_stats()     # environment readings from the instrument logs

mapper.save("prep_results/")   # both frames plus a record of how they were made
```

`save` writes `metadata.parquet`, `instrument_stats.parquet` and a
`provenance.jsonl` entry naming the package version, the configuration, which
logs were classified as what, and how many wells clashed. The metadata frame
on its own says none of that, and the classification in particular is worth
keeping: operators name labelling logs inconsistently, so it is the record of
what was actually picked up.

### AAS: substitution discovery

Every stage reads the same parameter file (see `examples/aas/params.yaml`) and
is run in order:

There is a manual database search in the middle: detection writes a FASTA
that has to be searched against the raw files before validation has anything
to read. So the pipeline runs in two phases:

```python
import proteolyzer.aas as aas

pipeline = aas.Pipeline("params.yaml")

pipeline.run_detection()    # preprocess -> translate -> detect
#   ... search the raw files against <output>/<sample>_validation.fasta ...
pipeline.run_validation()   # preprocess -> validate -> quantify

pipeline.status()           # what has run, and what can run now
```

`Pipeline` exists for the ordering: the preprocessor runs again in phase two
to convert the validation searches, the six-frame translation is skipped when
its frames are already on disk, and phase two refuses to start before the
searches exist. The stages can still be driven individually:

```python
aas.Preprocessor.MaxQuant(params).run()
aas.FrameTranslator(params).run()
aas.Detection(params).run()
aas.Validation(params).run()
aas.Quantification(params).run()
```

Stages exchange frames as parquet under the output folder (`.p` files from
earlier versions are still read), and each run appends its resolved
parameters, timestamp and package version to
`<output folder>/provenance.jsonl`.

To read a run back without knowing the stage-internal file names:

```python
results = aas.Results.from_params(params)   # or aas.Results("out/")

results.samples                             # what is in there
results.summary()                           # rows per sample per step; NA where
                                            # a step did not run, so reading down
                                            # a column shows where it stopped
results.combined("quantified")              # every sample in one frame, with a
                                            # Sample column
results.provenance()                        # what produced it
```

The steps are named for their result rather than the file that holds it:
`candidates`, `alt`, `filtered`, `fasta_entries`, `validated`, `evidence`,
`quantified`.

### Nomenclature

| term | meaning |
|---|---|
| **SAAP** | a peptide carrying an amino acid substitution |
| **BASE** | the unmodified peptide a SAAP is measured against |
| **ALT** | a mass shift with an alternative explanation: a known modification |

Outputs live under `SAAP/` and `ALT/`, and the quantification columns are
`SAAP.Sum`, `BASE.Sum`, `SAAP.Plex.<n>`, `BASE.Plex.<n>`. Folders written
before this naming (`MTP/`, `PTM/`) are still read back by `Results`. The
`PTM ppm` parameter is now `ALT ppm`; the old spelling still loads.

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
        models.py     Data, LoadedData, ProcessedData
        processor.py  dtypes, derived columns, labelling information
        matrix.py     pivot to a quantitative matrix
        io.py         parquet interchange between pipeline stages
        pipeline.py   stage plumbing: parameters, progress, provenance
    plots/            optional: plotting base class and relational plots
    cellenone/        optional: cellenONE export parsing and well mapping
    aas/              optional: amino acid substitution pipeline
    unimod/           optional: UniMod SQL queries, built on demand
tests/                pytest suite
examples/             runnable notebooks per module
docs/                 Docusaurus site; docs/docs is generated API reference
```

The core and `reference` are imported with the package. Everything else loads
on first access, so a missing extra only fails for the module that needs it —
`tests/test_package_boundaries.py` enforces that, including that importing
proteolyzer does not pull in matplotlib.

## Development

```bash
pip install -e '.[dev]'
pre-commit install                  # run the lint hooks on commit
pytest                              # test suite
ruff check .                        # lint
ruff format src tests               # format
mypy                                # type check (scope in pyproject.toml)
pydoc-markdown pydoc-markdown.yml   # regenerate docs/docs API reference
```

`make help` lists the same targets.

Tests that need an optional dependency skip themselves, so the suite runs on a
core install. Run it from the repository root.

## Contributing

Contributions are welcome! Please fork the repository, make your changes, and
submit a pull request. Ensure `ruff check .` and `pytest` pass, and add tests
for new behaviour.
