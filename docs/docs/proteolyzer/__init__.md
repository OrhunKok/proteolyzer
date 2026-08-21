---
sidebar_label: proteolyzer
title: proteolyzer
---

Proteolyzer: processing, analysis and visualization of proteomics data.

The core (``config``, ``utils``, ``transformers``) is imported eagerly. The
optional subpackages listed in ``__lazy__`` are imported on first attribute
access so that a missing extra (``pip install proteolyzer[aas]``) only fails
for the module that needs it.

## importlib

## config

## transformers

## utils

## Data

## configure\_logging

#### \_\_all\_\_

#### \_\_lazy\_\_

Subpackages imported on first attribute access, see :func:`__getattr__`.

#### \_\_getattr\_\_

```python
def __getattr__(name: str)
```

#### \_\_dir\_\_

```python
def __dir__() -> list[str]
```
