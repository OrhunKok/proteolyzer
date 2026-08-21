"""The entry point: turn a path into a :class:`Report`."""

from .models import Data, Report


def read(source, **kwargs) -> Report:
    """Read search-engine output into a :class:`Report`.

    ``source`` is a path, a string path, or a file-like object. Remaining
    keyword arguments configure :class:`Data`, e.g. ``load_all_columns=True``
    or ``extra_cols_to_load=["Q.Value"]``.

        report = proteolyzer.read("report.parquet")
        processed = report.process()
    """
    return Data(source=source, **kwargs).load()
