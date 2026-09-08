# Engineering notes

What was measured, what it cost, and what turned out not to be true. These are the
long form of things the code cannot carry without becoming unreadable: row counts
off a real export, a benchmark and its confounds, an assumption that shipped and
was wrong.

The division of labour, so nothing has to be written twice:

| where | what belongs there |
| --- | --- |
| a docstring | what the thing does, what it takes, what it hands back, and any gotcha that changes how you call it |
| [DECISIONS.md](https://github.com/OrhunKok/proteolyzer/blob/master/DECISIONS.md) | the claim, in a sentence or two — an index |
| these notes | the evidence, the numbers, the history, the things tried and dropped |
| [CHANGELOG.md](../changelog.md) | what a given version changed, and what a consumer has to do |

A docstring that starts explaining *why* for more than a line or two belongs here
instead, with a pointer left behind.

## The notes

- [Recognising a format](recognising-a-format.md) — how a file is identified, and
  the four column names that can never be part of it
- [Spectronaut](spectronaut.md) — everything measured off a real export, including
  the two things that shipped wrong first
- [Numbers and gaps](numbers-and-gaps.md) — what counts as missing, why a word
  does not, and where pandas loses digits
- [Performance](performance.md) — where the time in a read actually goes, and one
  optimisation that was measured and thrown away
