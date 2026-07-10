# MAGeCK2: Model-based Analysis of CRISPR screens

[![CI](https://github.com/davidliwei/mageck2/actions/workflows/ci.yml/badge.svg)](https://github.com/davidliwei/mageck2/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/mageck2.svg)](https://pypi.org/project/mageck2/)
[![License: BSD-3-Clause](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](LICENSE)

MAGeCK2 is a new computational tool, built upon MAGeCK, for CRISPR screening analysis.

Compared with MAGeCK, MAGeCK2 highlights the following additional functionalities:

* Paired sample analysis
* Processing and analyzing screens with UMIs
* Processing and analyzing paired-guide screens

See [ROADMAP.md](ROADMAP.md) for the planned direction of the project.

MAGeCK2 is compatible with the lagacy [MAGeCK](https://sourceforge.net/projects/mageck/). Everything is unchanged -- just need to switch the command from "mageck" to "mageck2".

# Installation

MAGeCK2 requires Python &ge; 3.7 and a C++ compiler (used to build the bundled RRA
component; g++ on Linux or clang++ on macOS, plus `make`). The Python dependencies
(numpy, scipy, pandas, matplotlib, statsmodels) are installed automatically.

Install the latest release from [PyPI](https://pypi.org/project/mageck2/):

    pip install mageck2

To install the latest development version directly from GitHub:

    pip install git+https://github.com/davidliwei/mageck2.git

or from a local clone:

    git clone https://github.com/davidliwei/mageck2.git
    cd mageck2
    pip install .

For detailed instructions, see [mageck2-doc](https://github.com/davidliwei/mageck2-doc).


# Quick start

Run a statistical test from a count table:

    mageck2 test -k count_table.txt -t treatment_samples -c control_samples -n output_prefix

Runnable example datasets are available in the [mageck2-demo](https://github.com/davidliwei/mageck2-demo) repository.


# Documentation

For instructions on installation, usage and running examples, see [mageck2-doc](https://github.com/davidliwei/mageck2-doc).


# MAGeCK

The older version of MAGeCK can be accessed via [sourceforge](https://sourceforge.net/projects/mageck/) or [bitbucket](https://bitbucket.org/liulab/mageck/src/master/).


# Contributing

Contributions are welcome. See [CONTRIBUTING.md](CONTRIBUTING.md) for how to get started,
[GOVERNANCE.md](GOVERNANCE.md) for how the project is run, and [ROADMAP.md](ROADMAP.md) for
where it is headed.


# Questions

Questions? Bug? Recommendations? Here are a few ways:

* Join our [Google group](https://groups.google.com/g/mageck?hl=en);
* Create an [issue](https://github.com/davidliwei/mageck2/issues) in the github repository. Issues from demo or documentations will need to be posted to [mageck2-demo](https://github.com/davidliwei/mageck2-demo) or [mageck2-doc](https://github.com/davidliwei/mageck2-doc), respectively.


# Version history

## 2026 — Version 0.2.0: Packaging modernization

* First release published to [PyPI](https://pypi.org/project/mageck2/): `pip install mageck2`.
* Migrated packaging to `pyproject.toml`; MAGeCK2 is now pip-installable on Python 3.7+ (including 3.12 and 3.13).
* Added continuous integration and a smoke-test suite.
* Added project governance and contribution documentation.

See [CHANGELOG.md](CHANGELOG.md) for the full list of changes.

## 2020.12.01 Version 0.1

* Add paired-guide support in count command
* Add UMI support in count command
* The source code released.
