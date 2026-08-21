"""UniMod export tests.

Building the database needs the network, so these tests work against a
hand-made SQLite file with the same shape as the real one.
"""

import sqlite3

import pandas as pd
import pytest

pytest.importorskip("sqlalchemy")
pytest.importorskip("xml2db")

from proteolyzer.unimod import UnimodDBLoader, UniModProcessor  # noqa: E402


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


def test_processor_exits_when_the_database_is_absent(tmp_path):
    with pytest.raises(SystemExit):
        UniModProcessor(db_file=str(tmp_path / "missing.db"))
