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

`make install` installs `.[dev]`, which is enough for `make lint` and
`make types`. `make test` also wants the `unimod` extra — without it,
`sqlalchemy` is missing and the unimod tests silently don't run — so run
`pip install -e '.[dev,unimod]'` first. `make docs` needs the `docs` extra,
`pip install -e '.[docs]'`, which `make install` does not pull in; CI's docs
job installs it separately for that reason, and a contributor running the
gate doesn't need mkdocs.

`make test`, `make lint` and `make types` are the gate: all three, green, before
a pull request. CI runs them on every push, the suite on three platforms, and
builds the docs `--strict` besides. `make test-downstream` is not part of it —
see below for why.

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
  making sure the consumers notice is not. Both of them run an `upstream.yml`
  workflow every Monday that compares their pin against this repository's latest
  release, moves it, runs their own suites against the new wheel, and opens a
  pull request saying whether they passed — a red pull request naming the
  release beats silence.
- **`make test-downstream`**, which runs their suites against this working tree
  when they are checked out under `downstream/`. That checkout is a session's,
  made by hand for the afternoon it's needed — not a CI fixture. **This
  repository is public and an Actions log is world-readable**, so running two
  private repositories' suites here would publish their test names, their paths
  and their failure output. The leak would not be a credential; it would be the
  console. That is why it is not wired in, and it does not stop being true when
  the credential problem is solved — the agent's App could clone them today.

  The integration signal exists already, from the side that can have it: each
  consumer's `upstream.yml` runs on Mondays, moves its pin to the newest release,
  runs its own suites against that wheel and opens a pull request saying whether
  they passed — in a private repository, where that output belongs.

  What this command adds over that is narrow and worth having: the working tree,
  before a release, rather than a published wheel after one. So it is a tool for
  a session to reach for when it is unsure, not a gate: it skips a consumer that
  isn't checked out rather than failing on it, which is right for the first and
  wrong for the second.

Breaking something on purpose is fine. Breaking it silently is what the tags are
against: bump the minor, write down what moved and what it is now, and say what a
consumer has to do.

## Reading another repository

There is no shared index, no mounted checkout and no session hook. If you need
to see how a consumer uses this package, clone it — it takes a second and needs
no setup:

```bash
gh repo clone OrhunKok/streamlit-DO-MS /tmp/x && grep -rn proteolyzer /tmp/x
```

`downstream/` is the same idea kept for an afternoon: check a consumer out there
and `make test-downstream` runs its suite against this working tree.

This repository is where the surviving copy of the cellenONE reader lives, and
that is worth knowing before writing a second one of anything: `CoordinatesMapping`
existed twice, in two repositories, and both copies fixed the same
imaging-channel bug in the same week without either knowing. Search before you
write — `gh api -X GET /search/code -f q="<term> user:OrhunKok"` — rather than
trusting an index to have been kept true.

## The agent that runs here

An issue labelled `agent` is worked by `.github/workflows/agent.yml`, in this
repository, using this repository's `CLAUDE.md` and `DECISIONS.md` as its brief.
Nothing outside this repository is involved, and nothing here reaches into
another repository: a change that belongs elsewhere is filed as an issue there
and worked by that repository's own agent.
