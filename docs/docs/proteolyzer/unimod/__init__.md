---
sidebar_label: unimod
title: proteolyzer.unimod
---

Query the UniMod database.

The database is not shipped: it is 5 MB against the 221 KB of CSVs the core
actually reads, so it is built once on first use and cached. Querying it needs
nothing beyond the standard library; *building* it downloads the UniMod XML and
needs the extra::

    pip install &#x27;proteolyzer[unimod]&#x27;

    from proteolyzer import unimod

    unimod.tables()                        # what is in there
    unimod.table(&quot;modifications&quot;)          # a whole table as a frame
    unimod.query(&#x27;&#x27;&#x27;
        SELECT m.full_name, s.one_letter, m.mono_mass
        FROM specificity AS s
        JOIN modifications AS m ON s.mod_key = m.record_id
        WHERE m.mono_mass BETWEEN ? AND ?
    &#x27;&#x27;&#x27;, (15.9, 16.1))

The cache lives under ``$PROTEOLYZER_CACHE_DIR``, else ``$XDG_CACHE_HOME``,
else ``~/.cache``. :func:`refresh` rebuilds it.

## os

## sqlite3

## Sequence

## Path

## pd

## PACKAGE\_LOGGER

#### \_\_all\_\_

#### DATABASE\_NAME

#### logger

#### cache\_dir

```python
def cache_dir() -> Path
```

Where the built database is kept.

#### database\_path

```python
def database_path() -> Path
```

The path the database is cached at, whether or not it exists yet.

#### database

```python
def database(build: bool = True) -> Path
```

The cached database, built on first use.

Parameters
----------
build
    Build it if it is missing. With ``False`` a missing database raises
    instead, which is what a caller that must not hit the network wants.

#### refresh

```python
def refresh(xml_source: str | None = None,
            xsd_source: str | None = None) -> Path
```

Download UniMod and rebuild the cached database.

Imported here rather than at module scope so that querying an existing
cache needs none of the build dependencies.

#### connect

```python
def connect(build: bool = True) -> sqlite3.Connection
```

A read-only connection to the cached database.

#### query

```python
def query(sql: str, params: Sequence = (), build: bool = True) -> pd.DataFrame
```

Run `sql` against the database and return the result as a frame.

#### tables

```python
def tables(build: bool = True) -> list[str]
```

The tables available to query.

#### table

```python
def table(name: str, build: bool = True) -> pd.DataFrame
```

A whole table as a frame.

