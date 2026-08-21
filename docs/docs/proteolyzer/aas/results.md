---
sidebar_label: results
title: proteolyzer.aas.results
---

Read back what the AAS pipeline produced.

The stages write one frame per sample per step, named after the step that
wrote them and split across ``SAAP/`` and ``ALT/``. That is fine for the
pipeline and unhelpful for the person asking &quot;what did we find?&quot;, who has to
know that the answer is ``SAAP/&lt;sample&gt;_SAAP_Quant.parquet`` and that stage 2
is validated while stage 1 is not.

Folders written before the SAAP/ALT/BASE naming (``MTP/``, ``PTM/``) are still
read: see :data:``0.

:class:``1 is the way in: it finds the samples, says which steps
completed for each, and combines a step&#x27;s frames across samples.

    from proteolyzer import aas

    results = aas.Results(&quot;out/&quot;)
    results.samples                     # what is in there
    results.summary()                   # rows per sample per step
    results.combined(&quot;quantified&quot;)      # every sample, with a Sample column
    results.provenance()                # how it was produced

## json

## Path

## pd

## frame\_exists

## read\_frame

## Logged

## PROVENANCE\_FILE

#### ARTEFACTS

What the stages write, under a name describing the result rather than the
step: (subdirectory, file stem). The order is the order of the pipeline.

#### LEGACY\_ARTEFACTS

Where the same artefacts lived before the SAAP/ALT/BASE naming, so an older
output folder can still be read back.

#### SAMPLE\_COL

Column the sample identifier is added under when frames are combined.

## Results Objects

```python
class Results(Logged)
```

The output folder of an AAS run.

#### \_\_init\_\_

```python
def __init__(output_dir: str | Path)
```

#### from\_params

```python
@classmethod
def from_params(cls, params) -> Results
```

Open the results of the run described by a parameter file or dict.

#### path

```python
def path(artefact: str, sample: str) -> Path
```

Where `artefact` for `sample` lives, whether or not it exists.

Falls back to the pre-SAAP/ALT layout if only that copy is there, so an
older results folder reads back unchanged.

#### \_path

```python
def _path(layout: dict, artefact: str, sample: str) -> Path
```

#### has

```python
def has(artefact: str, sample: str) -> bool
```

Whether the step that writes `artefact` got as far as `sample`.

#### samples

```python
@property
def samples() -> list[str]
```

Every sample with at least one artefact, in a stable order.

#### load

```python
def load(artefact: str, sample: str) -> pd.DataFrame
```

One sample&#x27;s frame for `artefact`.

#### combined

```python
def combined(artefact: str, samples: list[str] | None = None)
```

Every sample&#x27;s `artefact` in one frame, with a Sample column.

Samples the step did not reach are skipped and reported, so a partial
run still returns what it did produce.

#### summary

```python
def summary() -> pd.DataFrame
```

Rows per sample per step; NA where the step did not run.

Reading down a column shows where the pipeline stopped.

#### provenance

```python
def provenance() -> pd.DataFrame
```

What ran in this folder, one row per stage run, newest last.

