---
sidebar_label: formats
title: proteolyzer.core.formats
---

Input-format configuration for the core loading/processing pipeline.

Each search-engine block describes the files proteolyzer recognizes, the
columns worth loading from them, how those columns map onto the canonical
proteolyzer names, and which of them must stay numeric.

## dataclass

## field

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

#### MIN\_CATEGORICAL\_SAVING

Fraction of a column&#x27;s memory that turning it categorical has to save
for the conversion to be worth making. Set from measurement: on a real
report the columns that benefit save 49% or more and the ones that do
not save under 2% (or cost memory), so anything in between separates
them.

#### DIANN

#### MaxQuant

