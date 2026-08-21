---
sidebar_label: config
title: proteolyzer.utils.config
---

Input-format configuration for the core loading/processing pipeline.

Each search-engine block describes the files proteolyzer recognizes, the
columns worth loading from them, how those columns map onto the canonical
proteolyzer names, and which of them must stay numeric.

## dataclass

## field

## CONFIG

#### \_DIANN\_REPORT\_COLS

Columns loaded from the DIA-NN precursor reports. The main report and the
first-pass report share a schema, so they share this set.

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

#### EXCLUDE\_CAT\_CONVERSION

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

#### EXCLUDE\_CAT\_CONVERSION

## Config Objects

```python
@dataclass(frozen=True)
class Config()
```

#### COL\_MEDIAN\_THRESHOLD

#### CARDINALITY\_THRESHOLD

#### DIANN

#### MaxQuant

#### Protease

