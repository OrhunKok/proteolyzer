# proteolyzer

Processing, analysis and visualization of proteomics data. `proteolyzer.core`
recognizes search-engine output, reads it, narrows it and normalizes it onto one
schema; `reference` holds the domain constants; `cellenone`, `plots` and `unimod`
are optional subpackages imported on first use.

```bash
make install     # editable, with the dev extra
make test        # pytest
make lint        # ruff check + ruff format --check
make types       # mypy
make docs        # mkdocs build
```

`make test`, `make lint` and `make types` are the gate: all three, green, before
a pull request. CI runs them on every push, the suite on three platforms, and
builds the docs `--strict` besides.

## Why it is built this way

[DECISIONS.md](./DECISIONS.md) — the choices that are not visible from the code,
and the ones whose reasoning lives in a docstring worth finding.

## This is the repository others depend on

Two of them, and neither resolves it from PyPI:

| who | what they use |
|---|---|
| `streamlit-DO-MS` | `core.formats`, `core.loader`, `narrow`, `cellenone` — the whole reading and narrowing path |
| `decoder` | `reference.CODON_TABLE`, and `Data(...).load().frame` in one module |

**A change to a public interface reaches them before it reaches any test here.**
That happened on 2026-08-20: `utils` became `core` and the constants separated
from the settings, and decoder stopped importing at all for two days
(decoder#3). Hence:

- **Tags.** `CHANGELOG.md` says what each version broke, and
  `.github/workflows/release.yml` attaches a wheel to every tag. Consumers pin
  the wheel — a GitHub source archive carries no `.git`, so setuptools-scm falls
  back and every tag installs as `0.0.0`, which pip cannot compare.
- **`make test-downstream`**, which runs their suites against this working tree
  when they are checked out under `downstream/`: streamlit-DO-MS's four scripts,
  and decoder's `tests/test_imports.py`. A consumer this machine has no checkout
  of is skipped rather than failed, so the target is worth exactly as much as
  what is checked out beside it.

Breaking something on purpose is fine. Breaking it silently is what the tags are
against: bump the minor, write down what moved and what it is now, and say what a
consumer has to do.

## Ecosystem

This repository is one of several with a Claude devcontainer, which is what makes
it part of the network. The connective tissue is
[claude-shared](https://github.com/OrhunKok/claude-shared), in two places:
`/workspace-shared` is the read-only host mount, and it is only how the installer
gets found; `~/.claude-shared` is the working copy this container fast-forwards
by itself after every session. **Read the working copy** — the mount is a
checkout on the host and lags until a person pulls it, and a stale index is not
inert. The session hook names the path it read, prints who else is in the
network, and prints anything addressed here, at the start of every session.

**Before writing anything reusable, read `~/.claude-shared/CAPABILITIES.md`.**
It is keyed by capability, so it answers "has this already been solved?" The
other repositories are checked out read-only under `~/.claude-siblings`, so
reading one is `grep` rather than a question.

That index exists because of this repository: `CoordinatesMapping` was written
twice, in two places, and both copies fixed the same imaging-channel bug in the
same week without either knowing. The cellenONE reader here is the surviving one.

### Working across the repositories

The bus is GitHub issues on the *receiving* repository, labelled `from-upstream`
for something they need to know and `cross-repo-question` for something only they
can answer. Open means unhandled, closed is the acknowledgement.

**Filing one is the last resort, not the first.** A session here has push rights
on the other repositories and a clone is one command, so the default is to do the
work where it belongs — including moving a consumer's pin after a release here,
which is the half that is usually forgotten:

```bash
git clone https://github.com/OrhunKok/<other> /tmp/<other> && cd /tmp/<other>
# ... the change, in their style, with their tests ...
gh pr create --fill && gh pr merge --merge
```
