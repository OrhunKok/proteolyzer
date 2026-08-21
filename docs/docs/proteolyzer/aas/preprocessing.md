---
sidebar_label: preprocessing
title: proteolyzer.aas.preprocessing
---

Conversion of raw search-engine output into the parquet inputs of the pipeline.

Each search directory is filtered down to the columns and rows the downstream
stages need, then written next to the original text export as parquet.

## os

## Path

## np

## pd

## Stage

## Config

## column\_mapping

#### CONFIG

## Preprocessor Objects

```python
class Preprocessor()
```

Namespace grouping the per-search-engine preprocessors.

## MaxQuant Objects

```python
class MaxQuant(Stage)
```

#### FILES

#### FILES\_NEEDED

#### LOAD\_COLS

#### RENAME\_MAP

#### \_\_init\_\_

```python
def __init__(params, queue=None)
```

#### run

```python
def run()
```

#### \_search\_directories

```python
def _search_directories()
```

Every directory that may hold search output, `data_dir` included.

Yielding only the child names would skip a flat layout, where the
text files sit directly in the data folder.

#### \_convert\_directory

```python
def _convert_directory(subdir: Path) -> int
```

Converts the needed search files in `subdir`; returns how many.

#### \_filter

```python
def _filter(base: str, df: pd.DataFrame) -> pd.DataFrame
```

Apply the per-file row filtering and derived columns.

#### \_evidence

```python
def _evidence(df: pd.DataFrame) -> pd.DataFrame
```

#### \_allpeptides

```python
def _allpeptides(df: pd.DataFrame) -> pd.DataFrame
```

#### \_msms

```python
def _msms(df: pd.DataFrame) -> pd.DataFrame
```

#### \_peptides

```python
def _peptides(df: pd.DataFrame) -> pd.DataFrame
```

## DIANN Objects

```python
class DIANN()
```

#### FILES

#### FILE\_EXT

#### \_\_init\_\_

```python
def __init__(data_folder: str)
```

#### run

```python
def run()
```

