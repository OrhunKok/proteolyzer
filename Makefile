.PHONY: help install test test-downstream lint types format docs docs-serve clean

help:
	@echo "install  editable install with dev tooling"
	@echo "test     run the test suite"
	@echo "test-downstream  run streamlit-DO-MS and decoder's suites against this core"
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
# this suite -- it breaks one in theirs. Neither consumer's suite is pytest,
# and each has its own way of picking up the working tree, so each is wired in
# by name rather than discovered generically. Both are private repositories:
# a consumer this machine has no checkout of under downstream/ is skipped, not
# failed.
test-downstream:
	@status=0; \
	if [ -d downstream/streamlit-DO-MS ]; then \
		printf '\n=== downstream/streamlit-DO-MS ===\n'; \
		( cd downstream/streamlit-DO-MS \
		  && pip install -r requirements.txt \
		  && pip install -e "$(CURDIR)" \
		  && for t in test_file_load test_plot_inputs test_plots_render test_cellenone_upload; do \
		         python "tests/$$t.py" || exit 1; \
		     done \
		) || status=1; \
	else \
		echo "No downstream/streamlit-DO-MS; skipping."; \
	fi; \
	if [ -d downstream/decoder ]; then \
		printf '\n=== downstream/decoder ===\n'; \
		( cd downstream/decoder && pip install -e "$(CURDIR)" && python tests/test_imports.py ) || status=1; \
	else \
		echo "No downstream/decoder; skipping."; \
	fi; \
	exit $$status

lint:
	ruff check .
	ruff format --check src tests scripts

types:
	mypy

format:
	ruff check --fix src tests scripts
	ruff format src tests scripts

docs:
	mkdocs build --strict

docs-serve:
	mkdocs serve

clean:
	rm -rf build dist site .pytest_cache .ruff_cache .mypy_cache htmlcov .coverage
	find . -name '__pycache__' -type d -prune -exec rm -rf {} +
	find src -name '*.egg-info' -type d -prune -exec rm -rf {} +
