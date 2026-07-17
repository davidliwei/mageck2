# Changelog

All notable changes to MAGeCK2 are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Bundled the `mageckGSEA` C++ helper (vendored from legacy MAGeCK, BSD
  licensed). It is compiled and installed onto `PATH` alongside `RRA`, so
  `pathway --method gsea` — the default — now works from a pip/source install.
- Pathway smoke tests (`pathway --method gsea` and `--method rra`, end to end)
  and a small GMT test fixture.

### Fixed

- `pathway` with its default `--method gsea` no longer fails; the `mageckGSEA`
  helper it shells out to is now built and shipped (previously absent, causing
  an opaque `FileNotFoundError` on a missing intermediate file). See issue #10.
- External helper commands (`RRA`, `mageckGSEA`) now fail fast with a clear,
  actionable error when they exit non-zero or are missing from `PATH`, instead
  of being silently ignored and crashing later on a file the helper never wrote.
- `mageckGSEA` no longer maps pathway genes that are absent from the rank file
  to the top-ranked gene (via `map::operator[]` inserting a default index of 0),
  which silently inflated enrichment scores. Absent genes are now skipped, and a
  pathway with no overlapping genes scores zero instead of reading out of bounds.
- `mageckGSEA` guards against divide-by-zero / NaN in the enrichment-score
  computation for degenerate pathways (empty overlap, all-zero scores, or a set
  covering the entire ranked list).
- `mageckGSEA` no longer reports a degenerate pathway (one covering the entire
  ranked list, so no enrichment score can be computed) as maximally significant.
  Previously the permutation test resampled the same full set, found no
  exceedance of the zero score, and returned `p_permutation=0`; such pathways
  now short-circuit to a p-value of 1.0 so they cannot distort sorting or FDR.
- Corrected the `count` help text for the paired-guide options `--pg-start`,
  `--pg-end`, `--pg-start-2`, and `--pg-end-2`, which was accidentally copied from
  the UMI options (it described the position "of UMI" and referenced `--umi-start`
  / `--umi-end`). It now describes the second guide and the correct coordinate
  reference for each `--pairguide` mode. Also fixed the `--pg-min-read` help, which
  stated "Default 2" while the actual default is 3. See issue #22.
- `pathway --method gsea` no longer ranks the gene-ranking file's header row as a
  phantom gene. The default input is a `*.gene_summary.txt` whose first row is a
  header; `mageckGSEA` used to parse it as a gene named `id` with score 0,
  inflating the gene count and shifting the ranking, percentile, and permutation
  math for every pathway. `mageckGSEA` now accepts a `-H` / `--skip-header`
  switch, and `pathway --method gsea` passes it so the header is dropped.

## [0.2.0] - 2026-07-10

Modernization of packaging, testing, and project governance, and the first
release published to PyPI (`pip install mageck2`). MAGeCK2 is distributed as a
source distribution that compiles the bundled RRA C++ helper on install; a C++
compiler and `make` are required.

### Added

- Published to PyPI: `pip install mageck2`. Releases are built and uploaded by a
  GitHub Actions workflow using PyPI Trusted Publishing (OIDC, no stored tokens).
- `MANIFEST.in` so the source distribution bundles the RRA C++ sources, Makefile,
  and launch script needed to compile on install.
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
- False paired-guide warnings when `--pg-pair-only` files include a header row.
  `--pg-pair-only` files are now split on tab (or comma for `.csv`), matching the
  count-table readers, so sgRNA ids containing spaces are no longer mis-parsed.
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

[Unreleased]: https://github.com/davidliwei/mageck2/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/davidliwei/mageck2/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/davidliwei/mageck2/releases/tag/v0.1.0
