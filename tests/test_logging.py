"""Importing the package must not reconfigure logging for the host application."""

import io
import logging

import pytest

from proteolyzer.core.loader import DataLoader
from proteolyzer.core.logging import PACKAGE_LOGGER, MetaLogging, configure_logging
from proteolyzer.core.matrix import MatrixBuilder


@pytest.fixture(autouse=True)
def _restore_logging():
    """Put the package logger back the way import left it."""
    yield
    configure_logging()


def test_records_do_not_reach_the_root_logger():
    """Regression: the package used to call logging.basicConfig on import."""
    root = logging.getLogger()
    captured = logging.StreamHandler(io.StringIO())
    root.addHandler(captured)
    try:
        PACKAGE_LOGGER.info("hello")
    finally:
        root.removeHandler(captured)

    assert captured.stream.getvalue() == ""
    assert not PACKAGE_LOGGER.propagate


def test_configure_logging_directs_records_to_a_stream():
    stream = io.StringIO()
    configure_logging(stream=stream)

    PACKAGE_LOGGER.info("a message")

    assert "a message" in stream.getvalue()
    assert "INFO" in stream.getvalue()


def test_configure_logging_replaces_previous_handlers():
    first, second = io.StringIO(), io.StringIO()
    configure_logging(stream=first)
    configure_logging(stream=second)

    PACKAGE_LOGGER.info("only once")

    assert len(PACKAGE_LOGGER.handlers) == 1
    assert "only once" in second.getvalue()


def test_level_can_be_raised_to_silence_the_package():
    stream = io.StringIO()
    configure_logging(level=logging.ERROR, stream=stream)

    PACKAGE_LOGGER.info("not shown")
    PACKAGE_LOGGER.error("shown")

    assert "not shown" not in stream.getvalue()
    assert "shown" in stream.getvalue()


def test_class_loggers_are_namespaced_under_the_package():
    """Every class logger sits under `proteolyzer`, so one call re-levels them all."""
    assert DataLoader.logger.name == "proteolyzer.core.loader.DataLoader"
    assert DataLoader.logger.name.startswith(f"{PACKAGE_LOGGER.name}.")


def test_a_class_defined_outside_the_package_keeps_its_own_namespace():
    class Example(metaclass=MetaLogging):
        pass

    assert Example.logger.name == f"{__name__}.Example"


def test_memory_check_is_skipped_unless_debug_is_enabled():
    stream = io.StringIO()
    configure_logging(level=logging.INFO, stream=stream)

    # Not a DataFrame: reaching the body would log an error, so silence proves
    # the DEBUG gate short-circuited before doing any work.
    MatrixBuilder(None)._memory_check(object())
    assert stream.getvalue() == ""

    configure_logging(level=logging.DEBUG, stream=stream)
    MatrixBuilder(None)._memory_check(object())
    assert "does not have a memory_usage method" in stream.getvalue()
