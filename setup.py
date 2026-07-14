#!/usr/bin/env python
"""MAGeCK2 build script.

Project metadata is declared in ``pyproject.toml``. This file only customizes
the build to compile the bundled C++ helpers, which MAGeCK2 invokes as the
``RRA`` and ``mageckGSEA`` commands at runtime. The compiled binaries are
installed onto the environment's ``bin`` directory (alongside the ``mageck2``
entry script) so they are found on ``PATH``.
"""

import subprocess
import sys

from setuptools import setup
from setuptools.command.build import build

# (display name, source directory) for each bundled C++ helper. Each directory
# has a Makefile whose default target builds ``bin/<binary>``.
HELPERS = [("RRA", "rra"), ("mageckGSEA", "gsea")]


class BuildWithHelpers(build):
    """Compile the bundled C++ helpers before the normal build steps."""

    def run(self):
        for name, srcdir in HELPERS:
            try:
                returncode = subprocess.call(["make"], cwd=srcdir)
            except OSError as exc:
                sys.exit(
                    "CRITICAL: could not run `make` to build %s: %s. "
                    "A C++ compiler (g++/clang++) and make are required." % (name, exc)
                )
            if returncode != 0:
                sys.exit(
                    "CRITICAL: error compiling the %s source code. "
                    "Please check your C/C++ compilation environment." % name
                )
        super().run()


setup(
    cmdclass={"build": BuildWithHelpers},
    scripts=["bin/mageck2"],
    data_files=[("bin", ["rra/bin/RRA", "gsea/bin/mageckGSEA"])],
)
