"""Logging configuration utilities for proteolyzer.

Every class built with :class:`MetaLogging` gets a ``logger`` named after its
import path (e.g. ``proteolyzer.core.loader.DataLoader``), so the whole
package can be re-levelled or silenced at once::

    import logging

    logging.getLogger("proteolyzer").setLevel(logging.WARNING)

Handlers are attached to the ``proteolyzer`` logger, never to the root logger:
importing a library must not reconfigure logging for the host application.
:func:`configure_logging` changes level, format or destination, and
propagation stays off so records do not reach handlers installed elsewhere.
"""

import logging
import sys
from typing import TYPE_CHECKING, TextIO

import pandas as pd

#: Root logger of the package; all class loggers live underneath it.
PACKAGE_LOGGER = logging.getLogger(__package__.split(".")[0])

DEFAULT_FORMAT = "%(levelname)s - %(asctime)s - %(funcName)s - %(name)s - %(message)s"
DEFAULT_DATE_FORMAT = "%Y-%m-%d %H:%M:%S"


def configure_logging(
    level: int = logging.INFO,
    fmt: str = DEFAULT_FORMAT,
    datefmt: str = DEFAULT_DATE_FORMAT,
    stream: TextIO | None = None,
) -> logging.Logger:
    """Send proteolyzer's log records to `stream`, replacing any previous setup.

    Called once on import so that progress messages show up in a notebook
    without any setup. Call it again to change the level or destination, or
    ``logging.getLogger("proteolyzer").handlers.clear()`` to go quiet.
    """
    handler = logging.StreamHandler(stream if stream is not None else sys.stderr)
    handler.setFormatter(logging.Formatter(fmt=fmt, datefmt=datefmt))

    # Detach only; the stream belongs to whoever passed it in.
    for existing in list(PACKAGE_LOGGER.handlers):
        PACKAGE_LOGGER.removeHandler(existing)

    PACKAGE_LOGGER.addHandler(handler)
    PACKAGE_LOGGER.setLevel(level)
    # Records are handled here, so don't hand them to the root logger as well.
    PACKAGE_LOGGER.propagate = False

    return PACKAGE_LOGGER


configure_logging()


class MetaLogging(type):
    """Metaclass attaching a per-class logger and a memory-usage helper."""

    def __new__(cls, name, bases, dct):
        def _memory_check(self, data: pd.DataFrame) -> None:
            """
            Logs the memory usage of the provided data at DEBUG level.

            Self-gates on the logger's effective level: if DEBUG isn't enabled
            for this class's logger, returns immediately without computing
            memory usage. This matters because memory_usage(deep=True) walks
            every cell of every object column — expensive on large frames —
            so we don't want to pay that cost just to discard the result.

            Args:
                data: The pandas DataFrame to check memory usage for.
            """
            if not self.logger.isEnabledFor(logging.DEBUG):
                return
            try:
                summed_bytes = sum(data.memory_usage(deep=True))
                total_mib = round(summed_bytes / (1024**2), 1)
                self.logger.debug(f"Total memory usage of data: {total_mib} MiB")
            except AttributeError:
                self.logger.error(
                    "The provided data does not have a memory_usage method. "
                    "Make sure the data is a Pandas DataFrame"
                )
            except Exception as e:
                self.logger.error(
                    f"An unexpected error occurred during memory check: {e}"
                )

        dct["_memory_check"] = _memory_check
        log_class = super().__new__(cls, name, bases, dct)
        module = dct.get("__module__", "")
        log_class.logger = logging.getLogger(f"{module}.{name}" if module else name)

        return log_class


class Logged(metaclass=MetaLogging):
    """Base class for the classes that log.

    Subclassing this is equivalent to using :class:`MetaLogging` directly — the
    metaclass is inherited — but it also declares what the metaclass injects,
    which a type checker cannot infer. ``__slots__`` is empty so subclasses
    keep theirs.
    """

    __slots__ = ()

    if TYPE_CHECKING:
        logger: logging.Logger

        def _memory_check(self, data: pd.DataFrame) -> None: ...
