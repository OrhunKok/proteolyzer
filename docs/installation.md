# Installation

Python 3.14 or newer.

```bash
git clone https://github.com/OrhunKok/proteolyzer
cd proteolyzer

pip install .                 # core
pip install '.[aas]'          # amino acid substitution pipeline
pip install '.[unimod]'       # only to *build* the UniMod cache
pip install -e '.[dev]'       # editable, with test and lint tooling
pip install -e '.[docs]'      # this documentation
```

## Extras

| extra | what it adds |
|---|---|
| `aas` | the amino acid substitution pipeline |
| `unimod` | building the UniMod database; *querying* it needs nothing |
| `dev` | pytest, ruff, mypy, pre-commit |
| `docs` | mkdocs-material and mkdocstrings |

The core install pulls in pandas, numpy and pyarrow, and nothing else.

!!! note "quickdna"

    `aas` installs `quickdna` from a pinned wheel, because PyPI has no build
    for the supported interpreter. It is needed only by
    [`FrameTranslator`][proteolyzer.aas.translation.FrameTranslator]; the rest
    of the AAS pipeline imports and runs without it.

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
