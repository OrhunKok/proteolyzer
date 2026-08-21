---
sidebar_label: build
title: proteolyzer.unimod.build
---

## os

## sqlite3

## tempfile

## Any

## requests

## SmallInteger

## DataModel

## UnimodDBLoader Objects

```python
class UnimodDBLoader()
```

Builds the UniMod SQLite database from the published XML.

Only needed to (re)build the cache: querying an existing database needs
nothing beyond the standard library. See :func:`proteolyzer.unimod.refresh`.

#### XSD\_SOURCE

#### XML\_SOURCE

#### UNIMOD\_XML\_CONFIG

#### data\_model

#### \_\_init\_\_

```python
def __init__(db_output: str,
             xsd_source: str | None = None,
             xml_source: str | None = None) -> None
```

Initializes the DataModel, fetching XSD content and using a temp file.

#### \_fetch\_content

```python
def _fetch_content(source: str) -> str
```

Fetches content from a URL or reads from a local file path.

#### load\_and\_clean

```python
def load_and_clean() -> None
```

Executes the full process: fetching/parsing XML, insertion, and cleanup.

#### \_insert\_data

```python
def _insert_data() -> None
```

Fetches XML content and handles parsing/insertion using xml2db.

#### \_cleanup\_tables

```python
def _cleanup_tables() -> None
```

Drops &#x27;unimod*&#x27; tables and renames tables ending in &#x27;_row&#x27;.

