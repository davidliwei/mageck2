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
so the tests always describe this checkout.

Note this override applies in CI too: ``pip install .`` runs ``make`` in the
source tree, so the build outputs are present there alongside the copies
``data_files`` installs. The installed copies would then never be exercised, and
a regression in helper installation could pass CI even though an installed
MAGeCK2 cannot find its helpers. To keep that coverage, the PATH as it stood
before this hook ran is published as ``MAGECK2_INSTALLED_PATH``, and the tests
asserting the helpers are installed look there rather than at the build tree.
"""

import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

# (binary name, directory the Makefile writes it to)
HELPER_BIN_DIRS = [
    ("RRA", REPO_ROOT / "rra" / "bin"),
    ("mageckGSEA", REPO_ROOT / "gsea" / "bin"),
]


#: The PATH before this hook ran -- i.e. where an *installed* MAGeCK2 finds its
#: helpers. Published so tests can assert on the installation rather than on the
#: build tree that gets prepended below.
INSTALLED_PATH_ENV = "MAGECK2_INSTALLED_PATH"


def pytest_configure(config):
    """Prepend repo-built helper directories to PATH, most specific first."""
    os.environ.setdefault(INSTALLED_PATH_ENV, os.environ.get("PATH", ""))

    prefix = []
    for name, bindir in HELPER_BIN_DIRS:
        if (bindir / name).is_file():
            prefix.append(str(bindir))

    if prefix:
        os.environ["PATH"] = os.pathsep.join(prefix + [os.environ.get("PATH", "")])
