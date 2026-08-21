---
sidebar_label: translation
title: proteolyzer.aas.translation
---

Six-frame translation of a genome FASTA into searchable protein strings.

## pickle

## np

## SeqIO

## Stage

## FrameTranslator Objects

```python
class FrameTranslator(Stage)
```

Translates every FASTA entry in all six frames.

Two versions of each frame are written: the translation as-is, and one with
isoleucine rewritten to leucine, since the two are isobaric and cannot be
distinguished by mass.

#### ISOLEUCINE\_LEUCINE

str.translate table collapsing I onto L.

#### \_\_init\_\_

```python
def __init__(params, queue=None)
```

#### run

```python
def run()
```

#### \_count\_entries

```python
def _count_entries()
```

#### \_translate\_sequences

```python
def _translate_sequences()
```

#### \_write\_outputs

```python
def _write_outputs()
```

