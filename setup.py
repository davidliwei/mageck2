#!/usr/bin/env python
"""MAGeCK2 build script.

Project metadata is declared in ``pyproject.toml``. This file only customizes
the build to compile the bundled RRA C++ helper, which MAGeCK2 invokes as the
``RRA`` command at runtime. The compiled binary is installed onto the
environment's ``bin`` directory (alongside the ``mageck2`` entry script) so it
is found on ``PATH``.
"""

import subprocess
import sys

from setuptools import setup
from setuptools.command.build import build


class BuildWithRRA(build):
    """Compile the RRA C++ helper before the normal build steps."""

    def run(self):
        try:
            returncode = subprocess.call(["make"], cwd="rra")
        except OSError as exc:
            sys.exit(
                "CRITICAL: could not run `make` to build RRA: %s. "
                "A C++ compiler (g++/clang++) and make are required." % exc
            )
        if returncode != 0:
            sys.exit(
                "CRITICAL: error compiling the RRA source code. "
                "Please check your C/C++ compilation environment."
            )
        super().run()


setup(
    cmdclass={"build": BuildWithRRA},
    scripts=["bin/mageck2"],
    data_files=[("bin", ["rra/bin/RRA"])],
)
