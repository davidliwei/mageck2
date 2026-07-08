# Changelog

All notable changes to MAGeCK2 are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

Modernization of packaging, testing, and project governance. This work will be
tagged as the next release.

### Added

- `pyproject.toml` packaging (setuptools backend); MAGeCK2 is now pip-installable
  on Python 3.7+, including 3.12 and 3.13.
- Continuous integration (GitHub Actions) running the test suite across
  Python 3.9–3.13 on Linux and Python 3.13 on macOS.
- A `pytest` smoke-test suite verifying package import, RRA availability on
  `PATH`, the CLI, and an end-to-end `mageck2 test` run.
- Project governance and community documentation: `ROADMAP.md`, `GOVERNANCE.md`,
  `MAINTAINERS.md`, `CODEOWNERS`, `CONTRIBUTING.md`, and `CODE_OF_CONDUCT.md`.
- CI and license badges, pip-based installation, and a quick-start example in
  the README.

### Changed

- Migrated packaging from `setup.py` (which relied on the removed `distutils`
  and `build_py_2to3`) to `pyproject.toml`.
- The bundled RRA component now builds with the portable `c++` compiler driver
  (g++ on Linux, clang++ on macOS) instead of a hardcoded `g++`, so clang-only
  systems install cleanly.
- Declared previously-missing runtime dependencies explicitly: numpy, scipy,
  pandas, matplotlib, and statsmodels.
- Refreshed the README, including MAGeCK2 compatibility information.

### Fixed

- Installation failure on Python ≥ 3.12 caused by the removed `distutils` /
  `build_py_2to3` machinery.
- `NameError` when an explicit `--trim-5` value is supplied to the count command.
- `ValueError` / crash in `mle` on count tables containing spaces in sgRNA or
  gene names.
- Error when generating count reports in the RMD file.
- CNV correction failure with modern NumPy/Python when reading gene symbols from
  copy-number tables.
- Clear error for the disabled experimental `mle --bayes` option instead of an
  import crash or silent success.

### Removed

- Committed ctags artifacts (`mageck2/tags`, `rra/src/tags`); now gitignored.

## [0.1.0] - 2020-12-01

Initial public source release, built on
[MAGeCK](https://sourceforge.net/projects/mageck/).

### Added

- Paired-sample analysis.
- UMI-aware counting.
- Paired-guide (dual-sgRNA) screen counting.
- Core functionality carried forward from MAGeCK: sgRNA counting, quality
  control, RRA and MLE essentiality testing with copy-number correction, and
  functional enrichment.

[Unreleased]: https://github.com/davidliwei/mageck2/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/davidliwei/mageck2/releases/tag/v0.1.0
