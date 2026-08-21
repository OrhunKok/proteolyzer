"""Query the UniMod database.

The database is not shipped: it is 5 MB against the 221 KB of CSVs the core
actually reads, so it is built once on first use and cached. Querying it needs
nothing beyond the standard library; *building* it downloads the UniMod XML and
needs the extra::

    pip install 'proteolyzer[unimod]'

    from proteolyzer import unimod

    unimod.tables()                        # what is in there
    unimod.table("modifications")          # a whole table as a frame
    unimod.query('''
        SELECT m.full_name, s.one_letter, m.mono_mass
        FROM specificity AS s
        JOIN modifications AS m ON s.mod_key = m.record_id
        WHERE m.mono_mass BETWEEN ? AND ?
    ''', (15.9, 16.1))

The cache lives under ``$PROTEOLYZER_CACHE_DIR``, else ``$XDG_CACHE_HOME``,
else ``~/.cache``. :func:`refresh` rebuilds it.
"""

import os
import sqlite3
from collections.abc import Sequence
from pathlib import Path

import pandas as pd

from proteolyzer.core.logging import PACKAGE_LOGGER

__all__ = ["database", "database_path", "query", "refresh", "table", "tables"]

DATABASE_NAME = "unimod.db"

logger = PACKAGE_LOGGER.getChild("unimod")


def cache_dir() -> Path:
    """Where the built database is kept."""
    override = os.environ.get("PROTEOLYZER_CACHE_DIR")
    if override:
        return Path(override)
    xdg = os.environ.get("XDG_CACHE_HOME")
    base = Path(xdg) if xdg else Path.home() / ".cache"
    return base / "proteolyzer"


def database_path() -> Path:
    """The path the database is cached at, whether or not it exists yet."""
    return cache_dir() / DATABASE_NAME


def database(build: bool = True) -> Path:
    """The cached database, built on first use.

    Parameters
    ----------
    build
        Build it if it is missing. With ``False`` a missing database raises
        instead, which is what a caller that must not hit the network wants.
    """
    path = database_path()
    if path.exists():
        return path
    if not build:
        raise FileNotFoundError(
            f"No UniMod database at {path}. Call proteolyzer.unimod.refresh() "
            "to build it (needs the 'unimod' extra and network access)."
        )
    return refresh()


def refresh(xml_source: str | None = None, xsd_source: str | None = None) -> Path:
    """Download UniMod and rebuild the cached database.

    Imported here rather than at module scope so that querying an existing
    cache needs none of the build dependencies.
    """
    from .build import UnimodDBLoader

    path = database_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        path.unlink()

    logger.info(f"Building the UniMod database at {path}")
    UnimodDBLoader(
        db_output=str(path), xml_source=xml_source, xsd_source=xsd_source
    ).load_and_clean()
    return path


def connect(build: bool = True) -> sqlite3.Connection:
    """A read-only connection to the cached database."""
    path = database(build=build)
    return sqlite3.connect(f"file:{path}?mode=ro", uri=True)


def query(sql: str, params: Sequence = (), build: bool = True) -> pd.DataFrame:
    """Run `sql` against the database and return the result as a frame."""
    with connect(build=build) as connection:
        return pd.read_sql_query(sql, connection, params=tuple(params))


def tables(build: bool = True) -> list[str]:
    """The tables available to query."""
    with connect(build=build) as connection:
        rows = connection.execute(
            "SELECT name FROM sqlite_master WHERE type = 'table' ORDER BY name"
        )
        return [name for (name,) in rows]


def table(name: str, build: bool = True) -> pd.DataFrame:
    """A whole table as a frame."""
    available = tables(build=build)
    if name not in available:
        raise ValueError(f"No table '{name}'. Available: {available}")
    # Interpolated because a table name cannot be a bound parameter; the name
    # is checked against the catalogue above.
    return query(f'SELECT * FROM "{name}"', build=build)
