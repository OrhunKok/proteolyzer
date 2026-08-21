---
sidebar_label: validation
title: proteolyzer.aas.validation
---

Fragment-level validation of the substitutions proposed by detection.

## np

## pd

## frame\_exists

## read\_frame

## write\_frame

## Stage

## Validation Objects

```python
class Validation(Stage)
```

Keeps only candidates supported by fragment ions spanning the substitution.

#### \_\_init\_\_

```python
def __init__(params, queue=None)
```

#### process\_sample

```python
def process_sample(sample)
```

#### frags\_containing\_aas

```python
def frags_containing_aas(row: pd.Series)
```

Fragment indices whose ions span the substituted residue.

#### frag\_count

```python
def frag_count(row: pd.Series, frag_ev_merge: pd.DataFrame)
```

Count matched b/y fragments supporting the substitution in `row`.

Kept for use on a single candidate. :meth:`fragment_evidence` computes
the same numbers for a whole frame without rescanning the fragments
once per candidate.

#### \_FRAGMENT\_KEYS

Columns identifying one observed fragment.

#### fragment\_evidence

```python
def fragment_evidence(saap: pd.DataFrame,
                      frag_ev_merge: pd.DataFrame) -> pd.Series
```

Matched b/y fragment count per candidate, for the whole frame.

Equivalent to :meth:`frag_count` per row, but by joining the fragments
a candidate *would* produce against the ones observed. Filtering per
candidate instead costs one pass over the fragments each, which is
quadratic in the size of the run.

#### \_expected\_fragments

```python
def _expected_fragments(saap: pd.DataFrame) -> pd.DataFrame
```

One row per (candidate, scan, ion) the candidate should produce.

#### saap\_validate

```python
def saap_validate(val_evidence: pd.DataFrame, val_msms: pd.DataFrame,
                  saap: pd.DataFrame)
```

