"""Known idiosyncrasies of a search engine's own raw columns.

Each of these reads a column the way a specific engine actually writes it --
before, or entirely without, proteolyzer's rename onto the canonical schema.
streamlit-DO-MS loads with ``rename=False`` and its pages are written against
each engine's own column names, so a helper that only worked after renaming
would be unusable there.

Moved from streamlit-DO-MS, which had five variants of this spread over four
files -- including three separate reimplementations of pulling a run name out
of a path, none aware of the other two. See DECISIONS.md.
"""

from pathlib import PurePosixPath

import pandas as pd

#: The prefix FragPipe writes onto an entry from its contaminant library. A
#: contaminant is recognised by the name it was searched under, rather than a
#: hardcoded keratin list that goes stale as new instruments and reagents show
#: up in the data.
CONTAMINANT_MARK: str = "contam_"

#: What Philosopher wraps a ``Spectrum File`` entry in: the pepXML it wrote the
#: identification into, named after the mzML ``Spectrum`` was read from.
FRAGPIPE_FILE_PREFIX: str = "interact-"
FRAGPIPE_FILE_SUFFIX: str = ".pepXML"


def run_name_from_path(path: str, *, prefix: str = "", suffix: str = "") -> str:
    """The run name in a file path: its basename, minus a known decoration.

    Several engines record the run as a full path or a decorated file name
    rather than the bare name a report is grouped by, and stripping a
    directory off one and a prefix and suffix off another is the same
    operation once the basename is taken first. One helper, rather than a copy
    of it per engine.
    """
    name = PurePosixPath(str(path).replace("\\", "/")).name
    if suffix and name.endswith(suffix):
        name = name[: -len(suffix)]
    if prefix and name.startswith(prefix):
        name = name[len(prefix) :]
    return name


def fragpipe_run_names(spectrum: pd.Series, spectrum_file: pd.Series) -> pd.Series:
    """The run each FragPipe PSM came from.

    Philosopher writes ``Spectrum`` as ``(file).(scan).(scan).(charge)``, and a
    run name can itself contain dots, so the last three fields are taken off
    the end rather than splitting at the first one. A row with no ``Spectrum``
    falls back to ``Spectrum File``, stripped of the ``interact-`` prefix and
    ``.pepXML`` suffix Philosopher wraps it in.
    """
    # Filled rather than left null: an all-missing column has nothing for the
    # string ops below to infer a dtype from, and comes back floating-point
    # instead of the object series `.str[0]` expects. Overwritten by `.where`
    # regardless of what it becomes.
    from_spectrum = spectrum.fillna(".").astype(str).str.rsplit(".", n=3).str[0]
    from_file = spectrum_file.astype(str).apply(
        lambda path: run_name_from_path(
            path, prefix=FRAGPIPE_FILE_PREFIX, suffix=FRAGPIPE_FILE_SUFFIX
        )
    )
    return from_spectrum.where(spectrum.notna(), from_file)


def psm_identified(sequence: pd.Series) -> pd.Series:
    """Whether a MaxQuant scan carries an identified PSM.

    An unidentified scan's ``Sequence`` comes back as a single space, an empty
    string, or null, depending on the MaxQuant version and the parser that
    read it -- all three mean the same thing, and this is the one place that
    has to know it.
    """
    return sequence.fillna("").astype(str).str.strip() != ""


def decoy_flag(is_decoy: pd.Series) -> pd.Series:
    """JMod's ``is_decoy``, normalized to bool however it was written.

    A parquet round-trip keeps it a real bool; a csv one leaves the word
    ``True``/``False``, or the number ``1``/``0``, as text. A decoy is a decoy
    under any of them.
    """
    if pd.api.types.is_bool_dtype(is_decoy) or pd.api.types.is_numeric_dtype(is_decoy):
        return is_decoy.astype(bool)
    return is_decoy.astype(str).str.strip().str.casefold().isin({"true", "1"})


def contaminant_flag(protein: pd.Series, mark: str = CONTAMINANT_MARK) -> pd.Series:
    """Whether `protein` names an entry from FragPipe's contaminant library.

    FragPipe prefixes its contaminant library, so a contaminant is recognised
    by the name it was searched under rather than by a keratin list that goes
    stale.
    """
    return protein.astype(str).str.startswith(mark)
