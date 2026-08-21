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
- `count` writes `<prefix>.pg_unmapped.txt` when `--unmapped-to-file` is given in
  paired-guide mode, recording pairs that were counted but excluded from
  `pg_count.txt` — because the second guide matched no entry in `--list-seq-2`,
  or because the combination is absent from `--pg-pair-only`. Previously these
  were discarded with no record, so a filtered pair was indistinguishable from
  one that was never sequenced. The columns match `pg_count.txt`; a second guide
  that matched no library entry has no ID or gene, so it is identified by its
  sequence and its gene is reported as `None`. See issue #27.
- End-to-end tests for paired-guide counting, covering the allowed-pair filter,
  the new pair-level unmapped file, and read pairs whose first guide is unmatched.

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
- Library rows dropped for duplicating an earlier row's sequence are now reported
  actionably. The message names the library file, the dropped sgRNA IDs, and the
  representative ID that now carries their reads, instead of only counting them
  ("There are N sgRNAs with duplicated sequences"). It is also no longer emitted
  when nothing was dropped. See issue #27.
- The `--pg-pair-only` warning for an ID missing from a library now distinguishes
  the common cause from a genuine mismatch. When the ID was dropped as a duplicate
  sequence, the message says so, states that the pair will be excluded from the
  paired-guide counts, and names the representative ID to use instead. See issue #27.
- `--pg-pair-only` no longer silently excludes whole constructs from `pg_count.txt`
  when the pair file names an sgRNA ID that was dropped for duplicating an earlier
  sequence. Counting resolves a sequence to its one surviving ID, but the pair file
  was matched by the ID as written, so such pairs could never match even though
  their sequence pair had been counted correctly; users had to rewrite the pair
  file by hand. Those IDs are now resolved to the representative ID, with a log
  summary of how many entries were remapped and what they became. Where two pair-file
  entries resolve onto the same pair -- possible when both libraries contain
  duplicated sequences -- their reads are indistinguishable by sequence and are
  reported as one row, which is now warned about rather than silently summed.
  An sgRNA ID whose first row was dropped as a duplicate but which a later row
  re-admits under a sequence of its own is a real sgRNA, not an alias, and is no
  longer remapped. An ID that names more than one duplicated sequence has no
  single representative and is left unresolved with a warning, rather than
  resolved to whichever representative the library happened to list last.
  See issue #28.
- `count` no longer crashes with `TypeError` when writing the pair-level unmapped
  file for a run where some first reads matched no guide. In `--pairguide secondpair`
  mode the second-read window is recorded even when the first read matched nothing,
  under a `None` key; those entries are now skipped, and remain reported in the
  read-level `unmapped.txt`. See issue #27.

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
