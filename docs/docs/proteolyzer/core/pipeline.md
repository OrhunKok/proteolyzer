---
sidebar_label: pipeline
title: proteolyzer.core.pipeline
---

Shared plumbing for a pipeline stage.

A stage is one step of an analysis: it is constructed from resolved
parameters, reports progress on a queue, and records what it ran so a results
folder describes itself. None of that is specific to a particular pipeline, so
it lives here rather than in the subpackage that happened to need it first.

Subclasses add their own inputs and implement :meth:`run`.

## datetime

## json

## PackageNotFoundError

## version

## Path

#### PROVENANCE\_FILE

One JSON object per stage run, appended in the output folder.

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

Base class for a pipeline stage.

Parameters
----------
params
    The resolved parameters for the whole pipeline. A subclass that reads
    them from a file resolves them before calling this.
queue
    Optional queue receiving ``(stream, payload)`` progress messages.
    Defaults to a :class:`NullQueue`, which discards them.

#### \_\_init\_\_

```python
def __init__(params: dict, queue=None)
```

#### announce

```python
def announce() -> None
```

Report that the stage is ready, by name.

#### run

```python
def run() -> None
```

#### record\_run

```python
def record_run(output_dir: Path) -> Path
```

Append this stage&#x27;s parameters to the provenance log in `output_dir`.

Outputs are usually keyed only by sample, so a re-run with different
thresholds overwrites the previous results. This log is what makes an
output folder self-describing: which stage ran, when, from which
version, and with which resolved parameters.

#### package\_version

```python
def package_version() -> str
```

The installed proteolyzer version, or &quot;unknown&quot; from a bare checkout.

