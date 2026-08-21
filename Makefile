.PHONY: help install test test-downstream lint types format docs docs-serve clean

help:
	@echo "install  editable install with dev tooling"
	@echo "test     run the test suite"
	@echo "test-downstream  run the suites in downstream/ against this core"
	@echo "lint     ruff check + format check (what CI runs)"
	@echo "types    mypy over the modules listed in pyproject"
	@echo "format   apply ruff formatting and safe fixes"
	@echo "docs     build the documentation site (strict)"
	@echo "docs-serve  serve the documentation with live reload"
	@echo "clean    remove caches and build artefacts"

install:
	pip install -e '.[dev]'

test:
	pytest

# The pipelines that used to live here are their own repos now, so a change to
# core.io, core.logging, core.pipeline or reference no longer breaks a test in
# this suite -- it breaks one in theirs. Check those repos out under
# downstream/ and this runs them against the working tree.
test-downstream:
	@if [ ! -d downstream ]; then \
		echo "No downstream/ directory; nothing to check against."; exit 0; \
	fi
	@status=0; \
	for repo in downstream/*/; do \
		[ -d "$$repo/tests" ] || continue; \
		printf '\n=== %s ===\n' "$$repo"; \
		( cd "$$repo" && PYTHONPATH=src python -m pytest tests -q -p no:cacheprovider ) || status=1; \
	done; \
	exit $$status

lint:
	ruff check .
	ruff format --check src tests tools scripts

types:
	mypy

format:
	ruff check --fix src tests tools scripts
	ruff format src tests tools scripts

docs:
	mkdocs build --strict

docs-serve:
	mkdocs serve

clean:
	rm -rf build dist site .pytest_cache .ruff_cache .mypy_cache htmlcov .coverage
	find . -name '__pycache__' -type d -prune -exec rm -rf {} +
	find src -name '*.egg-info' -type d -prune -exec rm -rf {} +
