---
sidebar_label: params
title: proteolyzer.aas.params
---

Schema and loader for the AAS pipeline parameter file.

The YAML file is validated against :class:`ParamsSchema` and then flattened:
each section keeps its shared keys and is merged with the block named after
``Utils.Workflow``, so a stage can read ``params[&quot;Detection&quot;][&quot;ALT ppm&quot;]``
without knowing which search engine produced the data.

## Path

## Any

## Literal

## yaml

## BaseModel

## Field

## ParamsModel Objects

```python
class ParamsModel(BaseModel)
```

#### model\_config

## Utils Objects

```python
class Utils(ParamsModel)
```

#### Data\_Folder

#### Output\_Folder

#### Metadata\_File

#### Workflow

#### Labelling\_Setup

#### Label\_Plex

## Translation Objects

```python
class Translation(ParamsModel)
```

#### Genome\_FASTA

#### Translated\_Frames\_Folder

## DetectionMaxQuant Objects

```python
class DetectionMaxQuant(ParamsModel)
```

#### Detection\_PEP

#### AA\_Substitution\_ppm

#### ALT\_ppm

#### Positional\_Probability\_Threshold

#### C\_n\_term\_Modification\_Threshold

## Detection Objects

```python
class Detection(ParamsModel)
```

#### Protease

#### Protein\_FASTA

#### MaxQuant

## ValidationMaxQuant Objects

```python
class ValidationMaxQuant(ParamsModel)
```

#### Validation\_PEP

#### PIF

#### Fragment\_Evidence

## Validation Objects

```python
class Validation(ParamsModel)
```

#### MaxQuant

## QuantificationMaxQuant Objects

```python
class QuantificationMaxQuant(ParamsModel)
```

#### Minimum\_Quantity

## Quantification Objects

```python
class Quantification(ParamsModel)
```

#### MaxQuant

## ParamsSchema Objects

```python
class ParamsSchema(ParamsModel)
```

#### Utils

#### Translation

#### Detection

#### Validation

#### Quantification

#### \_load\_yaml

```python
def _load_yaml(filepath: str | Path) -> dict[str, Any]
```

Load raw parameters from a YAML file using safe_load.

#### load\_params

```python
def load_params(params: str | Path | dict) -> dict
```

Validate parameters and merge the workflow-specific sections.

Accepts a path to a YAML file or an already-loaded mapping.

