---
sidebar_label: config
title: proteolyzer.aas.config
---

Configuration for the AAS pipeline: input files, columns and label schemes.

## dataclass

## field

## ClassVar

## DIANN Objects

```python
@dataclass(frozen=True)
class DIANN()
```

#### FILES

#### LOAD\_COLS

Columns to load per file. ``None`` loads every column.

#### FILE\_EXTENSIONS

#### COLS\_RENAME\_MAPPING

#### FILES\_NEEDED

Files each stage needs to find in a search output directory.

## MaxQuant Objects

```python
@dataclass(frozen=True)
class MaxQuant()
```

#### FILES

#### LOAD\_COLS

Columns to load per file. ``None`` loads every column.

#### FILE\_EXTENSIONS

#### COLS\_RENAME\_MAPPING

#### FILES\_NEEDED

Files each stage needs to find in a search output directory.

## TMT Objects

```python
@dataclass(frozen=True)
class TMT()
```

TMT reference tables.

Only ``MQ_TMT_MAP`` is consumed at the moment (by ``Stage.metadata``);
the other two are kept as reference data for TMT-specific workflows.

#### TMT\_PLEX\_MAP

fmt: skip

#### MQ\_TMT\_MAP

fmt: skip

#### MASS\_SHIFT

## Config Objects

```python
@dataclass(frozen=True)
class Config()
```

#### DIANN

#### MaxQuant

#### TMT

