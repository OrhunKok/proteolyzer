---
sidebar_label: base
title: proteolyzer.aas.base
---

The AAS pipeline&#x27;s shared stage behaviour.

:class:`proteolyzer.core.pipeline.Stage` supplies what any pipeline needs —
parameters, a progress queue, provenance. This adds what is specific to this
one: parameters come from the AAS parameter file, and the samples to process
come from an experiment metadata spreadsheet.

## cached\_property

## Path

## np

## pd

## PROVENANCE\_FILE

## NullQueue

## \_CoreStage

## Config

## load\_params

#### CONFIG

#### \_\_all\_\_

#### METADATA\_COLS

Columns required from the experiment metadata spreadsheet.

## Stage Objects

```python
class Stage(_CoreStage)
```

Base class for the stages of the AAS pipeline.

Parameters
----------
params
    Path to a parameter YAML file, or an equivalent mapping.
queue
    Optional queue receiving ``(stream, payload)`` progress messages.

#### \_\_init\_\_

```python
def __init__(params, queue=None)
```

#### record\_run

```python
def record_run() -> Path
```

Log this run in the configured output folder.

#### metadata

```python
@cached_property
def metadata() -> pd.DataFrame
```

Experiment metadata, with TMT channels mapped to MaxQuant indices.

#### samples

```python
@property
def samples() -> np.ndarray
```

The unique sample IDs to process, in a stable order.

#### run

```python
def run() -> None
```

Process every sample in the metadata.

#### process\_sample

```python
def process_sample(sample: str) -> None
```

Handle one sample; implemented by each stage.

#### \_locate\_sample\_dir

```python
def _locate_sample_dir(sample: str, suffix: str = "") -> Path | None
```

Find the directory holding `sample`&#x27;s search output, if any.

