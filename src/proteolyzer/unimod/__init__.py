"""UniMod database retrieval and export.

Two steps, also exposed as the ``proteolyzer-unimod`` command:

    UnimodDBLoader   downloads the UniMod XML/XSD into a SQLite database
    UniModProcessor  exports amino acid and modification tables as CSV

The CSVs shipped in ``proteolyzer/resources`` were produced this way and back
the reference tables in :mod:`proteolyzer.config`.
"""

from .loader import UnimodDBLoader
from .processor import UniModProcessor

__all__ = ["UnimodDBLoader", "UniModProcessor"]
