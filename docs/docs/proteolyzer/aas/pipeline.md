---
sidebar_label: pipeline
title: proteolyzer.aas.pipeline
---

Run the AAS stages in order, in the two phases the method allows.

The five stages are not one automated run: detection writes a FASTA that has
to be searched against the raw files by the search engine, and only then do
the ``*_val`` directories the validation phase reads exist. So there are two
phases with a manual step between them:

    pipeline = aas.Pipeline(&quot;params.yaml&quot;)

    pipeline.run_detection()    # preprocess -&gt; translate -&gt; detect
    #   ... search the raw files against &lt;output&gt;/&lt;sample&gt;_validation.fasta ...
    pipeline.run_validation()   # preprocess -&gt; validate -&gt; quantify

    pipeline.status()           # what has run, what can run now

What this buys over calling the stages directly is the ordering: the
preprocessor runs again in the second phase to convert the validation
searches, translation is skipped when its frames already exist, and phase two
refuses to start before the searches are there.

## Path

## Logged

## Detection

## load\_params

## Preprocessor

## Quantification

## Results

## FrameTranslator

## Validation

#### VALIDATION\_SUFFIX

Suffix cellenONE-independent validation search directories carry.

## Pipeline Objects

```python
class Pipeline(Logged)
```

The AAS stages for one parameter file, run in the right order.

#### \_\_init\_\_

```python
def __init__(params, queue=None)
```

#### results

```python
@property
def results() -> Results
```

What this run has produced so far.

#### frames\_translated

```python
def frames_translated() -> bool
```

Whether the six-frame translation has already been written.

#### validation\_searches

```python
def validation_searches() -> list[Path]
```

The validation search directories, once the manual search has run.

#### status

```python
def status() -> dict
```

What has run and what can run now.

#### run\_detection

```python
def run_detection(translate: bool | None = None) -> Results
```

Phase one: preprocess, translate the genome, detect candidates.

Parameters
----------
translate
    Whether to run the six-frame translation. The default skips it when
    its frames are already on disk, since it is the slow step and its
    output only depends on the genome.

#### run\_validation

```python
def run_validation() -> Results
```

Phase two: convert the validation searches, validate, quantify.

