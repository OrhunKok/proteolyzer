---
sidebar_label: base
title: proteolyzer.aas.base
---

Shared plumbing for the AAS pipeline stages.

Every stage (preprocessing, translation, detection, validation,
quantification) is constructed the same way: from a parameter file or dict,
with an optional multiprocessing queue for progress reporting. :class:`Stage`
holds that plumbing so each stage module only contains its own analysis.

## datetime

## json

## cached\_property

## PackageNotFoundError

## version

## Path

## np

## pd

## Config

## load\_params

#### CONFIG

#### PROVENANCE\_FILE

One JSON object per stage run, appended in the output folder.

#### METADATA\_COLS

Columns required from the experiment metadata spreadsheet.

## NullQueue Objects

```python
class NullQueue()
```

A &#x27;do-nothing&#x27; queue to replace the real Queue when multiprocessing isn&#x27;t used.

#### put

```python
def put(item)
```

## Stage Objects

```python
class Stage()
```

Base class for the stages of the AAS pipeline.

Parameters
----------
params
    Path to a parameter YAML file, or an equivalent mapping.
queue
    Optional queue receiving ``(stream, payload)`` progress messages.
    Defaults to a :class:`NullQueue`, which discards them.

#### \_\_init\_\_

```python
def __init__(params, queue=None)
```

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

#### record\_run

```python
def record_run() -> Path
```

Append this stage&#x27;s parameters to the provenance log.

Outputs are keyed only by sample, so a re-run with different thresholds
overwrites the previous results. The log is what makes an output folder
self-describing: which stage ran, when, from which version, and with
which resolved parameters.

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

