"""Test configuration for MAGeCK2.

The suite shells out to the bundled C++ helpers by bare name (``RRA``,
``mageckGSEA``), so it exercises whichever build happens to come first on PATH.
That is the correct binary after ``pip install .`` -- setup.py compiles the
helpers and installs them alongside the package -- but not during ordinary local
development: an editable install never runs that build step, and a MAGeCK v1
package in the same environment ships commands with the *same names*. The suite
then silently tests v1's binaries and reports failures against fixes that are
present in this checkout, which is confusing and points nowhere useful.

Put this repository's own freshly built helpers first on PATH when they exist,
so the tests always describe this checkout. Falls back to the installed ones
when the helpers have not been built (e.g. in CI, where ``pip install .`` has
already placed the right binaries on PATH).
"""

import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

# (binary name, directory the Makefile writes it to)
HELPER_BIN_DIRS = [
    ("RRA", REPO_ROOT / "rra" / "bin"),
    ("mageckGSEA", REPO_ROOT / "gsea" / "bin"),
]


def pytest_configure(config):
    """Prepend repo-built helper directories to PATH, most specific first."""
    prefix = []
    for name, bindir in HELPER_BIN_DIRS:
        if (bindir / name).is_file():
            prefix.append(str(bindir))

    if prefix:
        os.environ["PATH"] = os.pathsep.join(prefix + [os.environ.get("PATH", "")])
