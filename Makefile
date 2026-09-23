.PHONY: help install test test-downstream lint types format docs docs-serve clean

# Tools come from .venv once `make install` has built one, and from PATH
# otherwise -- which is CI, whose workflows install into the runner's Python and
# call the tools directly. The path is spelled out rather than exported on PATH
# because the make macOS ships (3.81) looks a simple command up on the PATH it
# started with, so an exported .venv/bin is never searched.
BIN := $(if $(wildcard .venv/bin/python),$(CURDIR)/.venv/bin/)

help:
	@echo "install  build .venv with uv: editable, with the dev and unimod extras"
	@echo "test     run the test suite"
	@echo "test-downstream  run streamlit-DO-MS and decoder's suites against this core"
	@echo "lint     ruff check + format check (what CI runs)"
	@echo "types    mypy over the modules listed in pyproject"
	@echo "format   apply ruff formatting and safe fixes"
	@echo "docs     build the documentation site (strict)"
	@echo "docs-serve  serve the documentation with live reload"
	@echo "clean    remove caches and build artefacts"

install:
	uv venv --allow-existing
	uv pip install -e '.[dev,unimod]'

test:
	$(BIN)pytest

# A change to the reading path or to reference breaks a test in a consumer's
# suite before it breaks one here. Neither consumer's suite is pytest,
# and each has its own way of picking up the working tree, so each is wired in
# by name rather than discovered generically. Both are private repositories:
# a consumer this machine has no checkout of under downstream/ is skipped, not
# failed.
test-downstream:
	@status=0; \
	if [ -d downstream/streamlit-DO-MS ]; then \
		printf '\n=== downstream/streamlit-DO-MS ===\n'; \
		( cd downstream/streamlit-DO-MS \
		  && uv pip install --python "$(BIN)python" -r requirements.txt \
		  && uv pip install --python "$(BIN)python" -e "$(CURDIR)" \
		  && for t in test_file_load test_plot_inputs test_plots_render test_cellenone_upload; do \
		         "$(BIN)python" "tests/$$t.py" || exit 1; \
		     done \
		) || status=1; \
	else \
		echo "No downstream/streamlit-DO-MS; skipping."; \
	fi; \
	if [ -d downstream/decoder ]; then \
		printf '\n=== downstream/decoder ===\n'; \
		( cd downstream/decoder && uv pip install --python "$(BIN)python" -e "$(CURDIR)" && "$(BIN)python" tests/test_imports.py ) || status=1; \
	else \
		echo "No downstream/decoder; skipping."; \
	fi; \
	exit $$status

lint:
	$(BIN)ruff check .
	$(BIN)ruff format --check src tests scripts

types:
	$(BIN)mypy

format:
	$(BIN)ruff check --fix src tests scripts
	$(BIN)ruff format src tests scripts

docs:
	$(BIN)mkdocs build --strict

docs-serve:
	$(BIN)mkdocs serve

clean:
	rm -rf build dist site .pytest_cache .ruff_cache .mypy_cache htmlcov .coverage
	find . -name '__pycache__' -type d -prune -exec rm -rf {} +
	find src -name '*.egg-info' -type d -prune -exec rm -rf {} +
