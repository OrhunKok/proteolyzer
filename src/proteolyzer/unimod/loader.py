from xml2db import DataModel
from sqlalchemy import SmallInteger
from typing import Dict, Any, Optional
import sqlite3
import requests
import os
import tempfile
import argparse

class UnimodDBLoader:
    """
    A class to load UNIMOD XML data from a file path or URL into an SQLite database 
    using temporary files to handle content when needed.
    """
    
    XSD_SOURCE: str = "https://www.unimod.org/xmlns/schema/unimod_tables_1/unimod_tables_1.xsd"
    XML_SOURCE: str = 'https://www.unimod.org/xml/unimod_tables.xml'
    DB_OUTPUT: str = '../../../data/unimod.db' 
    DB_CONN: str = f"sqlite:///{DB_OUTPUT}"
    
    UNIMOD_XML_CONFIG: Dict[str, Any] = {
        "tables": {
            "unimod": {
                "fields": {
                    "majorVersion": {"type": SmallInteger}, 
                    "minorVersion": {"type": SmallInteger}
                }
            },
        },
    }

    data_model: DataModel
    
    def __init__(self, xsd_source: Optional[str] = None, xml_source: Optional[str] = None, db_conn: Optional[str] = None) -> None:
        """Initializes the DataModel, fetching XSD content and using a temp file."""
        
        self.xsd_source = xsd_source if xsd_source else self.XSD_SOURCE
        self.xml_source = xml_source if xml_source else self.XML_SOURCE
        connection_string: str = db_conn if db_conn else self.DB_CONN
        
        if db_conn and connection_string.startswith("sqlite:///"):
            self.DB_OUTPUT = connection_string.split("sqlite:///", 1)[1]
            
        xsd_content: str = self._fetch_content(self.xsd_source)
        
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.xsd', encoding='utf-8') as tmp_xsd:
            tmp_xsd.write(xsd_content)
            temp_xsd_path = tmp_xsd.name
            
        try:
            self.data_model = DataModel(
                xsd_file=temp_xsd_path,
                connection_string=connection_string,
                model_config=self.UNIMOD_XML_CONFIG
            )
        except Exception as e:
            print(f"Error initializing DataModel with XSD from {temp_xsd_path}: {e}")
            raise
        finally:
            os.remove(temp_xsd_path)
        
        print(f"Initialized DataModel for database: {connection_string}")
        
    def _fetch_content(self, source: str) -> str:
        """Fetches content from a URL or reads from a local file path."""
        if source.lower().startswith(('http://', 'https://')):
            print(f"Fetching content from URL: {source}...")
            try:
                response = requests.get(source, timeout=10)
                response.raise_for_status()
                return response.text
            except requests.exceptions.RequestException as e:
                raise IOError(f"Failed to fetch content from {source}: {e}")
        else:
            print(f"Reading content from local file: {source}...")
            with open(source, 'r', encoding='utf-8') as f:
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
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.xml', encoding='utf-8') as tmp_xml:
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
        conn: Optional[sqlite3.Connection] = None
        try:
            conn = sqlite3.connect(self.DB_OUTPUT)
            cursor: sqlite3.Cursor = conn.cursor()
            
            cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
            tables: list[tuple[str]] = cursor.fetchall() 
            
            changes_made: bool = False
            
            for table_tuple in tables:
                table_name: str = table_tuple[0]
                sql_cmd: Optional[str] = None
                
                if table_name.startswith('unimod'):
                    sql_cmd = f"DROP TABLE {table_name};"
                    
                elif table_name.endswith('_row'):
                    new_name: str = table_name.rstrip('_row')
                    sql_cmd = f"ALTER TABLE {table_name} RENAME TO {new_name};"

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


def main():
    """Parses command-line arguments and runs the UnimodDBLoader."""
    parser = argparse.ArgumentParser(
        description="Load UNIMOD XML data into an SQLite database."
    )
    
    parser.add_argument(
        '--xsd-source', 
        type=str, 
        default=UnimodDBLoader.XSD_SOURCE, 
        help=f"Source URL or path for the XSD schema file. Default: {UnimodDBLoader.XSD_SOURCE}"
    )
    parser.add_argument(
        '--xml-source', 
        type=str, 
        default=UnimodDBLoader.XML_SOURCE, 
        help=f"Source URL or path for the XML data file. Default: {UnimodDBLoader.XML_SOURCE}"
    )
    parser.add_argument(
        '--db-output', 
        type=str, 
        default=UnimodDBLoader.DB_OUTPUT, 
        help=f"Path to the UniMod SQLite database file. Default: {UnimodDBLoader.DB_OUTPUT}"
    )
    
    args = parser.parse_args()
    
    loader = UnimodDBLoader(
        xsd_source=args.xsd_source,
        xml_source=args.xml_source,
        db_conn=args.db_output
    )
    
    loader.load_and_clean()

if __name__ == "__main__":
    main() 