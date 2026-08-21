---
sidebar_label: proteolyzer
title: proteolyzer
---

Proteolyzer: processing, analysis and visualization of proteomics data.

The base suite is :mod:`proteolyzer.core` — recognizing search-engine output,
reading it, and normalizing it into a consistent frame — with
:mod:`proteolyzer.reference` holding the domain constants they share. Both are
imported eagerly.

Everything else is an optional subpackage, imported on first attribute access
so that a missing extra (``pip install proteolyzer[aas]``) only fails for the
module that needs it. ``tests/test_package_boundaries.py`` keeps that true.

## importlib

## core

## reference

## Data

## DataLoader

## DataProcessor

## MatrixBuilder

## configure\_logging

#### \_\_all\_\_

#### \_\_lazy\_\_

Subpackages imported on first attribute access, see :func:`__getattr__`.

#### \_MOVED

Modules that moved, kept importable so existing scripts keep working.

#### \_\_getattr\_\_

```python
def __getattr__(name: str)
```

#### \_\_dir\_\_

```python
def __dir__() -> list[str]
```

