.PHONY: help install test lint types format docs clean

help:
	@echo "install  editable install with dev tooling"
	@echo "test     run the test suite"
	@echo "lint     ruff check + format check (what CI runs)"
	@echo "types    mypy over the modules listed in pyproject"
	@echo "format   apply ruff formatting and safe fixes"
	@echo "docs     regenerate the API reference under docs/docs"
	@echo "clean    remove caches and build artefacts"

install:
	pip install -e '.[dev]'

test:
	pytest

lint:
	ruff check .
	ruff format --check src tests tools

types:
	mypy

format:
	ruff check --fix src tests tools
	ruff format src tests tools

docs:
	pydoc-markdown pydoc-markdown.yml

clean:
	rm -rf build dist .pytest_cache .ruff_cache .mypy_cache htmlcov .coverage
	find . -name '__pycache__' -type d -prune -exec rm -rf {} +
	find src -name '*.egg-info' -type d -prune -exec rm -rf {} +
