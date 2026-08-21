"""Tests for the UniMod query plugin.

Building the database needs the network, so these work against a hand-made
SQLite file with the same shape as the real one, pointed at through the cache
directory override.
"""

import sqlite3
from pathlib import Path

import pandas as pd
import pytest

pytest.importorskip("sqlalchemy")
pytest.importorskip("xml2db")

from proteolyzer import unimod  # noqa: E402
from proteolyzer.unimod.build import UnimodDBLoader  # noqa: E402
from proteolyzer.unimod.export import UniModProcessor  # noqa: E402


@pytest.fixture
def unimod_db(tmp_path):
    """A miniature UniMod database, as produced by xml2db before cleanup."""
    path = tmp_path / "unimod.db"
    with sqlite3.connect(path) as conn:
        conn.executescript(
            """
            CREATE TABLE amino_acids_row (
                one_letter TEXT, three_letter TEXT, full_name TEXT,
                num_H INT, num_O INT, num_C INT, num_N INT, num_S INT, num_Se INT
            );
            INSERT INTO amino_acids_row VALUES
                ('G', 'Gly', 'Glycine', 3, 1, 2, 1, 0, 0),
                ('L', 'Leu', 'Leucine', 11, 1, 6, 1, 0, 0),
                ('-', 'None', 'None', 0, 0, 0, 0, 0, 0);

            CREATE TABLE elements_row (element TEXT, full_name TEXT, mono_mass REAL);
            INSERT INTO elements_row VALUES
                ('H', 'Hydrogen', 1.007825035),
                ('O', 'Oxygen', 15.99491463),
                ('C', 'Carbon', 12.0),
                ('N', 'Nitrogen', 14.003074),
                ('S', 'Sulfur', 31.9720707),
                ('Se', 'Selenium', 79.9165196);

            CREATE TABLE modifications_row (
                record_id INT, mono_mass REAL, full_name TEXT, code_name TEXT,
                composition TEXT, username_of_poster TEXT, approved INT
            );
            INSERT INTO modifications_row VALUES
                (1, 15.994915, 'Oxidation', 'Oxidation', 'O', 'unimod', 1),
                (2, 42.010565, 'Acetylation', 'Acetyl', 'C2H2O', 'someone', 0);

            CREATE TABLE positions_row (record_id INT, position TEXT);
            INSERT INTO positions_row VALUES (1, 'Anywhere');

            CREATE TABLE classifications_row (record_id INT, classification TEXT);
            INSERT INTO classifications_row VALUES (1, 'Post-translational');

            CREATE TABLE specificity_row (
                one_letter TEXT, mod_key INT, position_key INT, classifications_key INT
            );
            INSERT INTO specificity_row VALUES ('G', 1, 1, 1), ('L', 2, 1, 1);

            CREATE TABLE unimod (majorVersion INT, minorVersion INT);
            INSERT INTO unimod VALUES (2, 8);
            """
        )
    return path


def _tables(path):
    with sqlite3.connect(path) as conn:
        return {
            row[0]
            for row in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")
        }


def test_cleanup_renames_row_tables_and_drops_unimod(unimod_db, monkeypatch):
    # Bypass __init__: it fetches the XSD over the network.
    loader = UnimodDBLoader.__new__(UnimodDBLoader)
    loader.db_output = str(unimod_db)

    loader._cleanup_tables()

    tables = _tables(unimod_db)
    assert "amino_acids" in tables
    assert "amino_acids_row" not in tables
    assert "unimod" not in tables
    # Regression: rstrip('_row') would also eat trailing characters of the name.
    assert "classifications" in tables
    assert "specificity" in tables


def test_processor_extracts_masses_and_adds_xle(unimod_db):
    UnimodDBLoader.__new__(UnimodDBLoader)
    loader = UnimodDBLoader.__new__(UnimodDBLoader)
    loader.db_output = str(unimod_db)
    loader._cleanup_tables()

    processor = UniModProcessor(db_file=str(unimod_db))
    processor.process_and_save()

    aas = pd.read_csv(processor.aa_output)
    mods = pd.read_csv(processor.mods_output)

    assert "-" not in aas["one_letter"].tolist()
    xle = aas.set_index("one_letter").loc["J"]
    assert xle["three_letter"] == "Xle"
    assert xle["mono_mass"] == pytest.approx(113.084064, abs=1e-5)
    # Only approved (or unimod-posted) modifications are exported.
    assert mods["full_name"].tolist() == ["Oxidation"]


def test_processor_reports_an_absent_database(tmp_path):
    """A library-style error, not sys.exit: the caller decides what to do."""
    with pytest.raises(FileNotFoundError, match="missing.db"):
        UniModProcessor(db_file=str(tmp_path / "missing.db"))


# ------------------------------------------------------------- the query API


@pytest.fixture
def cached_database(unimod_db, monkeypatch, tmp_path):
    """Point the cache at a prepared database, as if it had been built."""
    cache = tmp_path / "cache"
    cache.mkdir()
    monkeypatch.setenv("PROTEOLYZER_CACHE_DIR", str(cache))

    loader = UnimodDBLoader.__new__(UnimodDBLoader)
    loader.db_output = str(unimod_db)
    loader._cleanup_tables()

    (cache / unimod.DATABASE_NAME).write_bytes(unimod_db.read_bytes())
    return cache / unimod.DATABASE_NAME


def test_the_cache_location_is_overridable(monkeypatch, tmp_path):
    monkeypatch.setenv("PROTEOLYZER_CACHE_DIR", str(tmp_path))
    assert unimod.database_path() == tmp_path / "unimod.db"

    monkeypatch.delenv("PROTEOLYZER_CACHE_DIR")
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "xdg"))
    assert unimod.database_path().parent == tmp_path / "xdg" / "proteolyzer"


def test_tables_are_listed(cached_database):
    assert {"amino_acids", "modifications", "specificity"} <= set(unimod.tables())


def test_a_table_comes_back_as_a_frame(cached_database):
    aminos = unimod.table("amino_acids")
    assert set(aminos["one_letter"]) == {"G", "L", "-"}


def test_an_unknown_table_lists_the_real_ones(cached_database):
    with pytest.raises(ValueError, match="No table 'nope'"):
        unimod.table("nope")


def test_arbitrary_sql_with_parameters(cached_database):
    result = unimod.query(
        """
        SELECT m.full_name, s.one_letter
        FROM specificity AS s
        JOIN modifications AS m ON s.mod_key = m.record_id
        WHERE m.mono_mass BETWEEN ? AND ?
        """,
        (15.9, 16.1),
    )
    assert result["full_name"].tolist() == ["Oxidation"]
    assert result["one_letter"].tolist() == ["G"]


def test_the_database_is_opened_read_only(cached_database):
    with pytest.raises(sqlite3.OperationalError, match="readonly"):
        with unimod.connect() as connection:
            connection.execute("DROP TABLE amino_acids")


def test_a_missing_cache_is_reported_rather_than_built(monkeypatch, tmp_path):
    """`build=False` is for callers that must not reach the network."""
    monkeypatch.setenv("PROTEOLYZER_CACHE_DIR", str(tmp_path / "empty"))
    with pytest.raises(FileNotFoundError, match="refresh"):
        unimod.database(build=False)
    with pytest.raises(FileNotFoundError, match="refresh"):
        unimod.tables(build=False)


def test_refresh_delegates_to_the_builder(monkeypatch, tmp_path):
    """The build dependencies are imported only when a build happens."""
    monkeypatch.setenv("PROTEOLYZER_CACHE_DIR", str(tmp_path / "cache"))
    calls = {}

    class FakeLoader:
        def __init__(self, db_output, xml_source=None, xsd_source=None):
            calls["db_output"] = db_output
            calls["xml_source"] = xml_source

        def load_and_clean(self):
            Path(calls["db_output"]).write_bytes(b"")

    monkeypatch.setattr("proteolyzer.unimod.build.UnimodDBLoader", FakeLoader)

    path = unimod.refresh(xml_source="file:///somewhere.xml")

    assert path == unimod.database_path()
    assert calls["xml_source"] == "file:///somewhere.xml"
    assert path.exists()
