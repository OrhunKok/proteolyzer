"""Write the API reference into the build, one page per module.

Run by mkdocs-gen-files while the site builds, so the reference is never
committed and cannot drift from the source. The previous setup generated
markdown with pydoc-markdown and kept it in the repository, which needed a CI
job to catch the copy going stale.

mkdocstrings reads the source statically through griffe rather than importing
it, so the optional dependencies of the aas, cellenone and unimod plugins are
not needed to document them.
"""

from pathlib import Path

import mkdocs_gen_files

SOURCE = Path("src")
PACKAGE = SOURCE / "proteolyzer"
REFERENCE = Path("reference")

nav = mkdocs_gen_files.Nav()

for path in sorted(PACKAGE.rglob("*.py")):
    module_path = path.relative_to(SOURCE).with_suffix("")
    doc_path = path.relative_to(SOURCE).with_suffix(".md")
    parts = tuple(module_path.parts)

    if parts[-1] == "__init__":
        # A package's own docstring belongs on its section landing page, which
        # mkdocs-section-index folds into the section itself.
        parts = parts[:-1]
        doc_path = doc_path.with_name("index.md")
    elif parts[-1] == "__main__":
        continue

    nav[parts] = doc_path.as_posix()

    with mkdocs_gen_files.open(REFERENCE / doc_path, "w") as page:
        print(f"::: {'.'.join(parts)}", file=page)

    mkdocs_gen_files.set_edit_path(REFERENCE / doc_path, Path("..") / path)

with mkdocs_gen_files.open(REFERENCE / "SUMMARY.md", "w") as summary:
    summary.writelines(nav.build_literate_nav())
