---
sidebar_label: export
title: proteolyzer.unimod.export
---

## os

## sqlite3

## np

## pd

## UniModProcessor Objects

```python
class UniModProcessor()
```

Exports the reference CSVs that ship in ``proteolyzer/resources``.

A maintainer step: run it after :func:`proteolyzer.unimod.refresh` to pick
up a new UniMod release.

#### MODS\_OUTPUT

#### AA\_OUTPUT

#### \_\_init\_\_

```python
def __init__(db_file: str,
             mods_output: str | None = None,
             aa_output: str | None = None) -> None
```

Initializes the UniModDataProcessor with a database file path.

#### connect\_to\_db

```python
def connect_to_db() -> None
```

Establishes a connection to the SQLite database.

#### close\_connection

```python
def close_connection() -> None
```

Closes the database connection.

#### get\_modifications

```python
def get_modifications() -> pd.DataFrame
```

Fetches approved modification data.

#### get\_amino\_acids

```python
def get_amino_acids() -> pd.DataFrame
```

Fetches standard amino acid composition data.

#### get\_elements

```python
def get_elements() -> pd.DataFrame
```

Fetches element mono-mass data.

#### calculate\_aa\_masses

```python
def calculate_aa_masses(amino_acids_df: pd.DataFrame,
                        elements_df: pd.DataFrame) -> pd.DataFrame
```

Calculates mono-mass for amino acids and aggregates Isoleucine/Leucine into Xle (&#x27;J&#x27;).

#### process\_and\_save

```python
def process_and_save() -> None
```

Executes the entire data extraction, processing, and saving process.

