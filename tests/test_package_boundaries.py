"""The core must not drag the optional subpackages' dependencies in.

`import proteolyzer` has to work on a bare core install, and stay fast. These
run in a fresh interpreter because the test session itself has everything
imported already.
"""

import json
import subprocess
import sys
import textwrap

import pytest

#: Third-party modules only an optional subpackage needs.
OPTIONAL_DEPENDENCIES = (
    "Bio",
    "ahocorasick",
    "fastparquet",
    "openpyxl",
    "quickdna",
    "requests",
    "sqlalchemy",
    "xml2db",
    "yaml",
)

#: Core dependencies that are still only needed by one subpackage, so importing
#: the package should not pay for them.
DEFERRED_CORE_DEPENDENCIES = ("matplotlib", "seaborn", "scienceplots", "adjustText")


def _imported_after(statement: str, candidates) -> list[str]:
    """Which of `candidates` are in sys.modules after running `statement`."""
    code = textwrap.dedent(f"""
        import json, sys
        {statement}
        print(json.dumps([m for m in {list(candidates)!r} if m in sys.modules]))
    """)
    result = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=True
    )
    return json.loads(result.stdout)


def test_importing_the_package_pulls_no_optional_dependency():
    assert _imported_after("import proteolyzer", OPTIONAL_DEPENDENCIES) == []


def test_importing_the_package_does_not_pull_the_plotting_stack():
    """Plotting is lazy, so a script that never plots does not pay for it."""
    assert _imported_after("import proteolyzer", DEFERRED_CORE_DEPENDENCIES) == []


def test_the_core_pipeline_pulls_no_optional_dependency():
    statement = (
        "import proteolyzer as pz; "
        "d = pz.Data; m = pz.MatrixBuilder; r = pz.reference.amino_acid_masses()"
    )
    assert _imported_after(statement, OPTIONAL_DEPENDENCIES) == []


@pytest.mark.parametrize(
    ("subpackage", "expected"),
    [
        ("aas", "yaml"),
        ("cellenone", None),
        ("plots", "matplotlib"),
        # The build dependencies are imported by refresh(), not on import.
        ("unimod", None),
    ],
)
def test_subpackages_are_imported_on_first_access(subpackage, expected):
    """The other side of the boundary: touching one does load its stack."""
    candidates = (*OPTIONAL_DEPENDENCIES, *DEFERRED_CORE_DEPENDENCIES)
    loaded = _imported_after(f"import proteolyzer as pz; pz.{subpackage}", candidates)
    if expected is None:
        assert loaded == []
    else:
        assert expected in loaded


def test_lazy_attributes_are_discoverable():
    import proteolyzer as pz

    for name in ("aas", "cellenone", "plots", "unimod"):
        assert name in dir(pz)
        assert name in pz.__all__


def test_an_unknown_attribute_still_raises():
    import proteolyzer as pz

    missing = "nope"
    with pytest.raises(AttributeError, match="has no attribute 'nope'"):
        getattr(pz, missing)
