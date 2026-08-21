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
- **Plotting** — publication-styled relational plots (e.g. volcano plots).
- **cellenONE module** — maps single cells prepared on a cellenONE to well
  positions, and flags well/label clashes.
- **Alternate RNA decoding (AAS) module** — the pipeline used for discovery of
  amino acid substitutions and PTMs ([paper](https://decode.slavovlab.net/)).

## Installation

Python 3.14 or newer.

```bash
git clone https://github.com/OrhunKok/proteolyzer
cd proteolyzer
pip install .                 # core
pip install '.[aas]'          # amino acid substitution pipeline
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

report = pz.Data(source="report.parquet").load()   # DIA-NN report, recognized by name
processed = report.process()                       # dtypes, derived columns, labels

processed.unique_runs                              # {'run1', 'run2', ...}
processed.unique_ids                               # number of distinct precursors

matrix = (
    pz.MatrixBuilder(processed)
    .matrix_generation("Ms1.Area", index=["Precursor.Id"], columns=["Run"])
    .normalize_matrix(within_groups=["Run"], agg_func=np.nansum)
    .matrix
)
```

`Data` recognizes known DIA-NN/MaxQuant file names and reads only the columns
configured for them. Use `load_all_columns=True` for everything, or
`extra_cols_to_load=[...]` to add to the default set.

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
```

### AAS: substitution discovery

Every stage reads the same parameter file (see `examples/aas/params.yaml`) and
is run in order:

```python
import proteolyzer.aas as aas

params = "params.yaml"

aas.Preprocessor.MaxQuant(params).run()   # search output -> parquet
aas.FrameTranslator(params).run()         # six-frame genome translation
aas.Detection(params).run()               # candidate substitutions and PTMs
aas.Validation(params).run()              # fragment-level validation
aas.Quantification(params).run()          # substitution ratios
```

Stages exchange frames as parquet under the output folder (`.p` files from
earlier versions are still read), and each run appends its resolved parameters,
timestamp and package version to `<output folder>/provenance.jsonl`, so a
results folder records how it was produced.

### Reference data

Masses, the genetic code and protease rules live in one place, as immutable
mappings read on first use:

```python
from proteolyzer import reference

reference.amino_acid_masses()["A"]        # 71.03711...
reference.protease("Trypsin").allowed_counts
reference.CODON_TABLE["AUG"]              # "M"
```

They come from the UniMod CSVs in `proteolyzer/resources`. Regenerating those
is a maintainer job, so it lives outside the package and keeps its
dependencies out of the install:

```bash
pip install -r tools/unimod/requirements.txt
python -m tools.unimod load    --db-output tools/unimod/unimod.db
python -m tools.unimod process --db-file   tools/unimod/unimod.db \
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
tools/unimod/         maintainer script that regenerates resources/*.csv
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
