import os
import sqlite3
import tempfile
from typing import Any

import requests
from sqlalchemy import SmallInteger
from xml2db import DataModel


class UnimodDBLoader:
    """Builds the UniMod SQLite database from the published XML.

    Only needed to (re)build the cache: querying an existing database needs
    nothing beyond the standard library. See :func:`proteolyzer.unimod.refresh`.
    """

    XSD_SOURCE: str = (
        "https://www.unimod.org/xmlns/schema/unimod_tables_1/unimod_tables_1.xsd"
    )
    XML_SOURCE: str = "https://www.unimod.org/xml/unimod_tables.xml"
    UNIMOD_XML_CONFIG: dict[str, Any] = {
        "tables": {
            "unimod": {
                "fields": {
                    "majorVersion": {"type": SmallInteger},
                    "minorVersion": {"type": SmallInteger},
                }
            },
        },
    }
    data_model: DataModel

    def __init__(
        self,
        db_output: str,
        xsd_source: str | None = None,
        xml_source: str | None = None,
    ) -> None:
        """Initializes the DataModel, fetching XSD content and using a temp file."""

        self.db_output = db_output
        self.xsd_source = xsd_source if xsd_source else self.XSD_SOURCE
        self.xml_source = xml_source if xml_source else self.XML_SOURCE
        connection_string: str = f"sqlite:///{db_output}"

        xsd_content: str = self._fetch_content(self.xsd_source)

        with tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".xsd", encoding="utf-8"
        ) as tmp_xsd:
            tmp_xsd.write(xsd_content)
            temp_xsd_path = tmp_xsd.name

        try:
            self.data_model = DataModel(
                xsd_file=temp_xsd_path,
                connection_string=connection_string,
                model_config=self.UNIMOD_XML_CONFIG,
            )
        except Exception as e:
            print(f"Error initializing DataModel with XSD from {temp_xsd_path}: {e}")
            raise
        finally:
            os.remove(temp_xsd_path)

        print(f"Initialized DataModel for database: {connection_string}")

    def _fetch_content(self, source: str) -> str:
        """Fetches content from a URL or reads from a local file path."""
        if source.lower().startswith(("http://", "https://")):
            print(f"Fetching content from URL: {source}...")
            try:
                response = requests.get(source, timeout=10)
                response.raise_for_status()
                return response.text
            except requests.exceptions.RequestException as e:
                raise OSError(f"Failed to fetch content from {source}: {e}") from e
        else:
            print(f"Reading content from local file: {source}...")
            with open(source, encoding="utf-8") as f:
                return f.read()

    def load_and_clean(self) -> None:
        """Executes the full process: fetching/parsing XML, insertion, and cleanup."""
        self._insert_data()
        self._cleanup_tables()
        print("\nUNIMOD data loading and cleanup complete.")

    def _insert_data(self) -> None:
        """Fetches XML content and handles parsing/insertion using xml2db."""

        xml_content: str = self._fetch_content(self.xml_source)

        print("Inserting data into target tables...")
        with tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".xml", encoding="utf-8"
        ) as tmp_xml:
            tmp_xml.write(xml_content)
            temp_xml_path = tmp_xml.name

        try:
            document = self.data_model.parse_xml(xml_file=temp_xml_path)
            document.insert_into_target_tables()
            print("Data insertion finished.")
        except Exception as e:
            print(f"Error during XML parsing or data insertion: {e}")
            raise
        finally:
            os.remove(temp_xml_path)

    def _cleanup_tables(self) -> None:
        """Drops 'unimod*' tables and renames tables ending in '_row'."""
        print("\nStarting table cleanup...")
        conn: sqlite3.Connection | None = None
        try:
            conn = sqlite3.connect(self.db_output)
            cursor: sqlite3.Cursor = conn.cursor()

            cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
            tables: list[tuple[str]] = cursor.fetchall()

            changes_made: bool = False

            for table_tuple in tables:
                table_name: str = table_tuple[0]
                sql_cmd: str | None = None

                if table_name.startswith("unimod"):
                    sql_cmd = f'DROP TABLE "{table_name}";'

                elif table_name.endswith("_row"):
                    # removesuffix, not rstrip: rstrip strips any trailing
                    # '_', 'r', 'o' or 'w' characters, not the suffix.
                    new_name: str = table_name.removesuffix("_row")
                    sql_cmd = f'ALTER TABLE "{table_name}" RENAME TO "{new_name}";'

                if sql_cmd:
                    try:
                        cursor.execute(sql_cmd)
                        changes_made = True
                    except sqlite3.Error as e:
                        print(f"  Error executing SQL on {table_name}: {e}")

            if changes_made:
                conn.commit()
                print("Cleanup changes committed.")
            else:
                print("No cleanup changes required.")

        except sqlite3.Error as e:
            print(f"Error connecting to or operating on SQLite database: {e}")
        finally:
            if conn:
                conn.close()
