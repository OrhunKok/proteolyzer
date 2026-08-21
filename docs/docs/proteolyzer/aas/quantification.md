---
sidebar_label: quantification
title: proteolyzer.aas.quantification
---

Quantification of validated substitutions relative to their base peptide.

## np

## pd

## frame\_exists

## read\_frame

## write\_frame

## Stage

## Quantification Objects

```python
class Quantification(Stage)
```

Computes mistranslated-to-base-peptide intensity ratios per sample.

#### \_\_init\_\_

```python
def __init__(params, queue=None)
```

#### process\_sample

```python
def process_sample(sample)
```

#### \_apply\_minimum\_quantity

```python
def _apply_minimum_quantity(quant: pd.DataFrame) -> pd.DataFrame
```

Drop ratios where either peptide is below ``Minimum Quantity``.

#### \_raas

```python
def _raas(mtp, mtp_df, bp_df, sample_df, label_designation)
```

Relative abundance of the substituted peptide vs its base peptide.

