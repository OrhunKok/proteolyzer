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
  back and every tag installs as `0.0.0`, which pip cannot compare. What the
  release still owes by hand is the tag, the wheel, and that changelog entry;
  making sure the consumers notice is not. Both of them run a `pin-check`
  workflow that compares their pin against this repository's latest release
  every Monday and hands their own worker a pull request carrying the wheel URL
  and the notes.
- **`make test-downstream`**, which runs their suites against this working tree
  when they are checked out under `downstream/`. Neither is wired in yet; both
  have a suite worth running.

Breaking something on purpose is fine. Breaking it silently is what the tags are
against: bump the minor, write down what moved and what it is now, and say what a
consumer has to do.

## Ecosystem

This repository is one of several with a Claude devcontainer, which is what makes
it part of the network. The connective tissue is
[claude-shared](https://github.com/OrhunKok/claude-shared), mounted read-only at
`/workspace-shared`. The session hook prints who else is in it, and anything
addressed here, at the start of every session.

**Before writing anything reusable, read `/workspace-shared/CAPABILITIES.md`.**
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
work where it belongs. Moving a consumer's pin after a release here is no longer
that work — the `pin-check` workflow does it on its own — but a from-upstream
issue that needs code on their side still does:

```bash
git clone https://github.com/OrhunKok/<other> /tmp/<other> && cd /tmp/<other>
# ... the change, in their style, with their tests ...
gh pr create --fill && gh pr merge --merge
```
