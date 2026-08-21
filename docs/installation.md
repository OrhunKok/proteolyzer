# Installation

Python 3.14 or newer.

```bash
git clone https://github.com/OrhunKok/proteolyzer
cd proteolyzer

pip install .                 # core
pip install '.[unimod]'       # only to *build* the UniMod cache
pip install -e '.[dev]'       # editable, with test and lint tooling
pip install -e '.[docs]'      # this documentation
```

## Extras

| extra | what it adds |
|---|---|
| `unimod` | building the UniMod database; *querying* it needs nothing |
| `dev` | pytest, ruff, mypy, pre-commit |
| `docs` | mkdocs-material and mkdocstrings |

The core install pulls in pandas, numpy and pyarrow, plus the plotting stack.

## Domain pipelines

The cellenONE and amino acid substitution pipelines were moved into their own
repositories, so neither their dependencies nor their release cadence land on
everyone installing the core:

- [proteolyzer-cellenone](https://github.com/OrhunKok/proteolyzer-cellenone)
- [proteolyzer-aas](https://github.com/OrhunKok/proteolyzer-aas)

Each depends on this package and uses `core.io`, `core.logging` and
`core.pipeline`.

!!! note "LaTeX for the default plot theme"

    The plotting module defaults to the `science`
    ([scienceplots](https://github.com/garrettj403/SciencePlots)) theme, which
    renders text with LaTeX. Without a LaTeX installation, pass another theme:
    `VolcanoPlot(..., theme="default")`.

## Development

```bash
pip install -e '.[dev]'
pre-commit install       # run the lint hooks on commit

pytest                   # test suite
ruff check .             # lint
ruff format src tests    # format
mypy                     # type check (scope in pyproject.toml)
mkdocs serve             # documentation, with live reload
```

`make help` lists the same targets.

Tests that need an optional dependency skip themselves, so the suite runs on a
core install. Run it from the repository root.
