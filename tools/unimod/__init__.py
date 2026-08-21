"""Regenerate the UniMod tables that ship in ``proteolyzer/resources``.

A maintainer tool, not part of the library: it downloads the UniMod XML into a
SQLite database and exports the two CSVs the package reads. Users analysing
data never run it, which is why it lives outside the importable package and
keeps its dependencies out of the install.

    pip install -r tools/unimod/requirements.txt
    python -m tools.unimod load    --db-output tools/unimod/unimod.db
    python -m tools.unimod process --db-file   tools/unimod/unimod.db \
        --mods-output src/proteolyzer/resources/unimod_modifications.csv \
        --aa-output   src/proteolyzer/resources/unimod_amino_acids.csv
"""

from .loader import UnimodDBLoader
from .processor import UniModProcessor

__all__ = ["UnimodDBLoader", "UniModProcessor"]
