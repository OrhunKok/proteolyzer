---
sidebar_label: logging
title: proteolyzer.core.logging
---

Logging configuration utilities for proteolyzer.

Every class built with :class:`MetaLogging` gets a ``logger`` named after its
import path (e.g. ``proteolyzer.core.loader.DataLoader``), so the whole
package can be re-levelled or silenced at once::

    import logging

    logging.getLogger(&quot;proteolyzer&quot;).setLevel(logging.WARNING)

Handlers are attached to the ``proteolyzer`` logger, never to the root logger:
importing a library must not reconfigure logging for the host application.
:func:`configure_logging` changes level, format or destination, and
propagation stays off so records do not reach handlers installed elsewhere.

## logging

## sys

## TYPE\_CHECKING

## TextIO

## pd

#### PACKAGE\_LOGGER

Root logger of the package; all class loggers live underneath it.

#### DEFAULT\_FORMAT

#### DEFAULT\_DATE\_FORMAT

#### configure\_logging

```python
def configure_logging(level: int = logging.INFO,
                      fmt: str = DEFAULT_FORMAT,
                      datefmt: str = DEFAULT_DATE_FORMAT,
                      stream: TextIO | None = None) -> logging.Logger
```

Send proteolyzer&#x27;s log records to `stream`, replacing any previous setup.

Called once on import so that progress messages show up in a notebook
without any setup. Call it again to change the level or destination, or
``logging.getLogger(&quot;proteolyzer&quot;).handlers.clear()`` to go quiet.

## MetaLogging Objects

```python
class MetaLogging(type)
```

Metaclass attaching a per-class logger and a memory-usage helper.

#### \_\_new\_\_

```python
def __new__(cls, name, bases, dct)
```

## Logged Objects

```python
class Logged(metaclass=MetaLogging)
```

Base class for the classes that log.

Subclassing this is equivalent to using :class:`MetaLogging` directly — the
metaclass is inherited — but it also declares what the metaclass injects,
which a type checker cannot infer. ``__slots__`` is empty so subclasses
keep theirs.

#### \_\_slots\_\_

