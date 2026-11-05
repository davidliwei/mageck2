# Changelog

All notable changes to MAGeCK2 are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `RRA --max-permutation-per-group <n>`: an optional upper bound on the
  permutations performed for each distinct sgRNA-per-gene count, trading run
  time for p value resolution. There is no bound by default, and the option is
  ignored when `--permutation` is given explicitly. A 20,000-gene screen
  performs 2.02 million permutations per group and spends 2.7 s per RRA call
  doing it, twice per comparison in `mageck2 test`; bounding it at 500,000
  brings that to 0.76 s.

  The bound is opt-in because it is not free, and the cost is easy to miss.
  `calculateFDRVector` multiplies the smallest p value by the number of tested
  genes, and the permutation count scales with that same number, so for a gene
  that beats the entire null distribution the two cancel and its adjusted p
  value is `1 / (2 * (rand_passnum + 1))` — about 0.005 whatever the library
  size. A fixed bound breaks the cancellation: the smallest reachable adjusted
  p value becomes `(number of genes) / (2 * bound)`, which for 60,000 genes
  against a 500,000 bound is 0.06, so no gene in such a screen could be called
  at FDR 0.05 however strong it was. Screens with many hits do not show this,
  because the block of genes tied at the p value floor is pulled down by the
  step-down pass; a screen with one real hit shows it at once. Set the option
  with a target floor in mind — `(number of genes) / (2 * floor)` permutations
  keep an isolated hit callable at it.

- `tests/test_rra.py` additionally covers RRA's memory and name handling: peak
  RSS must scale with the screen rather than with a compiled-in worst case, a
  200-character name must survive intact, an over-long one must be truncated
  and reported rather than overflowing, a screen past the initial group
  capacity must come back complete, and an over-long name must still be matched
  by `--control` and by `--skip-gene`.

- `tests/test_rra.py`: tests for the bundled RRA helper's own contract. The
  suite reached RRA only through `mageck2 test` and `mageck2 pathway`, and
  asserted little beyond a zero exit status, so nothing pinned down how many
  permutations RRA performs, what p value resolution that buys, or what that
  resolution costs in achievable significance. The last of those is now an
  explicit invariant: the adjusted p value floor for an isolated hit is fixed
  at 1/202 and does not move with the library size.

- `mageckGSEA --seed <int>`: the seed for the permutation test, between 1 and
  2147483646, default 123456. The seed is reported in the run log, so an output
  file can be traced back to the stream that produced it. Varying it across runs
  is how to measure the Monte Carlo noise the permutation p values carry; the
  enrichment score and the analytic p value do not depend on it.

- `tests/test_gsea.py`: tests for the bundled mageckGSEA helper's permutation
  contract. The suite checked that scores and p values were well formed but
  never that a run was repeatable, which it was not. The regression test lets
  the clock tick over between the two runs it compares, because the old seeding
  read the clock in whole seconds and two runs inside the same second agreed
  even with the bug present.

- `tests/parity/` and `tests/test_rra_parity.py`: a differential-test harness
  for the planned replacement of the bundled C++ helpers with Python, together
  with goldens recorded from the current compiled RRA. The goldens can only be
  produced while the C++ exists, so they are recorded now; afterwards they
  remain as the oracle, and the last commit containing the C++ remains in git
  history as the way to regenerate them.

  Thirteen cases cover the paths a port is most likely to get wrong, including
  the two that no previous testing touched: the by-group permutation path,
  which is the default, and `--control`, which replaces the uniform null with
  resampled control-sgRNA percentiles. Two further cases round values coarsely
  so that sgRNAs tie by the dozen, because percentiles come from
  `bTreeSearchingF`, whose tie resolution follows the recursion path rather
  than any stated rule.

  The harness proves its own oracle and its own comparator rather than assuming
  either: goldens are written only when two consecutive runs of the helper
  agree, and the comparator is mutation-tested against known-wrong candidates,
  since a differential harness that cannot fail converts "we did not check"
  into "we checked and it passed". One case is deliberately left without a
  golden — `--control` with a list matching nothing reaches an out-of-bounds
  read (issue #51), and freezing a golden there would promote one build's
  undefined behaviour to a requirement.

- `tests/parity/gsea/` and `tests/test_gsea_parity.py`: the same harness for
  `mageckGSEA`, with goldens recorded from the current compiled helper. Twenty
  cases cover both of `parserankfile`'s unrelated scoring branches — `-c 0`
  derives a gene's score from its line number and had no coverage at all —
  along with the tie regimes that decide output ordering, and the four
  degenerate gene sets real GMT files contain: one covering the whole ranked
  list, one overlapping it not at all, one of duplicates, and one of a single
  gene.

  Recording the goldens turned up five input-handling bugs in the helper, each
  confirmed against the binary and each now a corpus case that requires the
  future port to differ rather than to agree:

  - A GMT file whose last line has no trailing newline loses its last gene set
    silently. `parsegmt` breaks out of its read loop on the same `getline` that
    returned the line.
  - The same in `parserankfile`, and worse: dropping the last gene also lowers
    the line count every gene's score is divided by, so one missing newline
    perturbs the whole ranking rather than losing one row.
  - A blank line anywhere in a GMT file is stored as a gene set named with the
    empty string. `main` then matches it against the default `--pathway_name`
    of `""`, takes the single-pathway branch, and reports that one empty set
    while silently skipping every real one — exiting 0.
  - `parsegmt` tokenises each line on whitespace rather than on tabs, although
    GMT is tab-delimited and its second column is free-text description. The
    description's second and later words are therefore read as gene names.
    Usually they match nothing and only inflate the reported `Size`, but gene
    symbols appear in prose descriptions routinely, and one that matches is
    added to the set: giving a two-gene set the description `response to G01`
    moved its enrichment score from 5.55e-17 to 0.396226 and its p value from
    0.9025 to 0.
  - A GMT line ending in whitespace gains a phantom empty gene name, because
    `while(ss.good())` runs one iteration past the last token and a failed
    extraction into a `std::string` leaves it empty rather than untouched. It
    inflates the reported `Size` by one.

  The enrichment score is not required to agree to the last bit, and the
  reason is worth recording: `getgseascore_fast` computes
  `decvalue*(ni-i) + sumincval`, which clang contracts into a fused
  multiply-add at `-O3` wherever FMA is in the baseline instruction set — it is
  on arm64, it is not on baseline x86-64. Rebuilding the same source with
  `-ffp-contract=off` changes three of the twenty goldens, by 2.776e-17 and
  5.551e-17, in gene sets whose score cancels to zero. So the helper does not
  agree with itself across compiler flags there, and requiring byte equality
  would have failed the C++ against its own goldens the first time CI ran on a
  different architecture. The permutation p values are still required to agree
  exactly: the 32-bit truncation of `gs_ret::gscore` absorbs the difference
  before the permutation comparison sees it.

- `mageck2.rra`: a Python implementation of RRA, with a documented array-based
  API — sgRNA names, gene names and values in, a result object out — rather
  than a file path in and a file path out. It is not yet the default: the
  compiled helper still is, and the Python path is selected by setting
  `MAGECK2_RRA_IMPL=python`. That switch is an internal migration control, not
  a user-facing option; it is read from the environment rather than added to
  the command line because offering the choice would invite users to make it,
  and there is no basis on which they could.

  The two implementations agree exactly, not merely closely. Across the
  thirteen recorded cases `lo_value`, `p` and `FDR` each deviate by zero at the
  precision the output file carries, with no change in gene ranking: the
  bundled Lehmer generator is reproduced bit for bit, so both implementations
  draw the same permutation null rather than two samples from the same
  distribution, and `scipy.special.betainc` agrees with the bundled
  `BetaNoncentralCdf` to about 1e-11. A 20,000-gene screen agrees equally
  exactly, and `mageck2 test` produces a byte-identical `gene_summary.txt`
  either way. The full test suite passes under both.

  On 20,000 genes × 5 sgRNAs with 2.02 million permutations the Python
  implementation uses 0.51 s of CPU against the helper's 0.44 s, and 240 MB of
  peak memory against 98 MB — 56 MB of which is the cost of importing
  NumPy/SciPy at all, which MAGeCK2 already pays. The permutation null is built
  in chunks for that reason; computing it in one array cost 541 MB on the same
  screen.

  Output differs from the helper's in one respect, which the comparison treats
  as immaterial: rows that tie on `lo_value` may appear in a different order,
  because the helper sorts them with an unstable quicksort and this sorts them
  stably. No value differs.

- `tests/test_rra_python.py`: tests for the Python implementation's own
  contract — the array API, the environment switch, the generator's behaviour
  across the block boundaries it draws in, and the refusal to run when a
  `--control` list matches nothing. The last of these is a chosen behaviour
  rather than a copied one: the helper reads out of bounds in that case, and
  falling back to a uniform null would report p values that are not the
  control-based ones the caller asked for while saying nothing about it.

### Changed

- RRA allocates for the screen it was given rather than for the largest one
  imaginable. `MAX_NAME_LEN` was 10,000 and is held inline by every
  `ITEM_STRUCT` and `GROUP_STRUCT`, so each gene cost 10,056 bytes and each
  sgRNA 10,040 whatever its name actually was. Peak RSS measured on a
  single-list screen:

  | genes | before | after |
  |---|---|---|
  | 2,000 | 123 MB | 11 MB |
  | 8,000 | 485 MB | 36 MB |
  | 20,000 | 1,209 MB | 87 MB |
  | 60,000 | 3,463 MB | 258 MB |

  `MAX_NAME_LEN` is now 512, which clears the longest names RRA is given
  (MSigDB gene set ids run to roughly 200 characters), and the three copies
  into those fields — unbounded `strcpy` before, safe only because the field
  was 10,000 bytes — are bounded, with over-long names truncated in the output
  and reported on stderr. Grouping is unaffected, since it keys on the full
  name; `--control` and `--skip-gene` are matched against the stored, truncated
  form, so they keep working at any name length. Two names that agree over the
  first 511 characters are indistinguishable to all three of those, which is
  what the stderr warning says. Output is byte-identical on every fixture
  tested, including ragged
  sgRNA counts, multiple lists, and comma-joined multi-gene rows. A 20,000-gene
  screen also runs about 25% faster (0.44 s to 0.33 s), from the smaller
  `memcpy`s in `QuickSortGroupByLoValue`.

  The groups array additionally starts at 4,096 genes and doubles as the input
  is counted, instead of reserving `MAX_GROUP_NUM` (100,000) entries up front.
  That reservation was lazily mapped, so it cost address space rather than
  resident memory, but it also meant the existing growth path was dead code —
  it did `delete[]` on a `malloc`-ed pointer and left the caller to `free` a
  `new[]`-ed one. It is now `realloc` throughout, and it is the normal path for
  any genome-scale screen.

- `argsParser.build_parser()` builds the command-line parser and returns it;
  `crisprseq_parseargs()` now calls it and then parses. No behavior change — the
  split exists so the parser can be inspected without side effects, which lets a
  documentation build generate the option reference from the parser itself rather
  than from `--help` output pasted in by hand. The hand-pasting is why `USAGE.md`
  still described `--pairguide auto` after 0.3.0 removed it.

- Reworded a few help strings to the phrasing the documentation had been using,
  where nothing was lost by it: `mle --threads`, and
  `mle --sgrna-eff-name-column`/`--sgrna-eff-score-column`, which now say the
  column index is 0-based rather than leaving it to be inferred from the default.

### Fixed

- `mageckGSEA` reads its two input files correctly. Five input-handling bugs,
  all of which produced plausible-looking output and none of which reported
  anything. The most consequential affects essentially every real GMT file:

  - **A gene set's size was overstated, and sometimes its score was wrong.**
    GMT is tab-delimited and its second column is a free-text description, but
    `parsegmt` split lines on *whitespace* and treated every token after the
    second as a gene. The description's second and later words therefore
    joined the gene set. Usually they match nothing and only inflate the
    `Size` column — the bundled `tests/data/pathways.gmt` has five genes per
    set and was reported as 7 — but gene symbols appear in prose descriptions
    constantly, and one that matches is *scored*: giving a two-gene set the
    description `response to G01` moved its enrichment score from 5.55e-17 to
    0.396226 and its p value from 0.9025 to 0, because `G01` was the
    top-ranked gene. Lines are now split on tabs.
  - **A blank line anywhere in a GMT file suppressed the entire analysis.**
    It was stored as a gene set named with the empty string, which then
    matched the default `--pathway-name` of `""`, so `main` took its
    single-set branch, reported that one empty set, and silently skipped every
    real set in the file — exiting 0, having written an output file that looks
    like a result. Blank lines are now ignored.
  - **A file not ending in a newline lost its last line.** Both readers used
    `while(true){ getline(...); if(eof()) break; }`, and the `getline` that
    consumes a final unterminated line also sets `eofbit`, so the loop broke
    before using what it had just read. In a GMT file that silently dropped
    the last gene set. In a ranked gene list it was worse than losing one row:
    the dropped gene also lowered the line count every gene's weight is
    divided by, so one missing newline perturbed the whole ranking.
  - **A GMT line ending in whitespace gained a phantom member**, because
    `while(ss.good())` runs one iteration past the last token and a failed
    extraction into a `std::string` leaves it empty rather than untouched.
    It matched nothing but inflated `Size` by one.

  Nothing changes for input that was already well formed. That is asserted
  rather than claimed: all fifteen well-formed cases in `tests/parity/gsea/`
  still reproduce their goldens — recorded from the *previous* build —
  bit for bit, and exactly the five malformed cases changed.

  The ranked list's weights are still divided by one more than the number of
  genes. The old loop ran once past the end and incremented that counter
  anyway, so correcting the loop would have changed the divisor too, and with
  it every number MAGeCK2 has ever produced. It is preserved deliberately;
  whether that normalisation should change is a separate decision from fixing
  the input handling, and is not made here.

  The five cases were previously classified in the parity harness as
  divergences — behaviour the Python port was required to *differ* from.
  Fixing the helper is the better answer while the helper is still what ships,
  so they are ordinary parity cases now, and the property they carry is
  required of both implementations: the answer on a malformed file must equal
  the answer on its repair. `mageckGSEA` and `mageck2.gsea` now produce
  byte-identical output on every case, the bundled GMT included.

- The Python RRA backend reads the group column as the list it is. RRA's second
  field is comma-separated and a row belongs to *every* name in it
  (`fileio.cpp:107`, `:266`), which is exactly how `mageck2 pathway --method
  rra` works: `pathwayFunc.py` writes one row per gene whose group field joins
  every pathway containing that gene. The backend read the field as a single
  name, so a gene in three pathways produced a pathway called `P1,P2,P3` that
  does not exist while the three real ones each lost a member. On a 6-pathway
  GMT with 62 of 100 genes in more than one set, the summary came back listing
  invented pathways such as `OVERLAP00,OVERLAP01` and omitting every real one.
  A repeated name is not collapsed either, matching the helper: `P1,P1` gives
  P1 the row twice.

  The bundled `pathways.gmt` has no gene in two pathways, so no test in the
  suite ever put a comma in that column. `tests/parity/` gains
  `pathway_memberships` and `pathway_duplicate_membership`, recorded from the
  helper on inputs in the shape the pathway caller really emits, and
  `tests/test_smoke.py` gains an end-to-end check over an overlapping GMT that
  asserts the identity of the reported pathways under both backends.

- The Python RRA backend groups on full names. It truncated group identifiers
  to `MAX_NAME_LEN - 1` *before* grouping, so two names agreeing to that length
  merged into one gene — changing its score and the size of the population the
  permutation and FDR are computed over. The helper keys `groupNames` on the
  untruncated string and truncates only the label it stores, which is what its
  own warning promises. Measured on a probe input: the helper prints two rows
  under one 511-character label with 4 and 6 sgRNAs; the backend printed one
  row of 10.

- The Python RRA backend no longer ignores the optional `prob` and `chosen`
  columns. Both change the answer in the helper — `prob` below 1 switches to
  `ComputeLoValue_Prob`, and `chosen` of 0 drops a row from its pool's
  percentiles and its gene's lo-value while still counting it in that gene's
  sgRNA total — and neither is implemented, so a file carrying a meaningful
  value was silently analysed as though the column were absent. Such a file is
  now refused with a message naming `MAGECK2_RRA_IMPL=cpp`. The six-column
  files `mageck2` writes itself always carry the inert values (`prob` is the
  literal 1 and only chosen rows are written), so no shipped path is affected;
  a test pins that.

- The Python RRA backend's permutation chunking is bounded by matrix elements
  rather than by rows, so a wide gene cannot defeat the memory limit. At the
  previous 2**18-row limit a gene with 100 sgRNAs materialized a 262,144 x 100
  float64 array — 210 MB, and the reduction holds several arrays of that shape
  at once — which is the input the limit exists for. Tiling and saturation
  libraries do put that many sgRNAs on one gene.

- `mageckGSEA` permutation p values are reproducible. The permutation test was
  seeded with `srand(time(NULL))` and drew with `rand()`, so two identical
  invocations more than a second apart disagreed: on a 50-pathway GMT against a
  5,000-gene ranking, `p_permutation` differed for 45 of the 50 pathways, `FDR`
  for all 50, and the reported ranking for 31. A pathway result could not be
  reproduced, including by whoever produced it. Permutations now run on the
  bundled Lehmer generator in `rngs.cpp` -- the one RRA already seeds with the
  same constant -- so a run is a function of its inputs and of `--seed`.
  `ES`, the analytic `p`, `Hits` and `LFC` are unaffected, since none of them
  were ever drawn from the generator; `p_permutation`, `FDR` and `Ranking` will
  not match what earlier versions reported, which is the point.

- `mageckGSEA` no longer reads one element past the end of the score vector
  while permuting. `rand()` may return `RAND_MAX`, and the index was computed as
  `int(rand()*1.0/RAND_MAX*nsz)`, which is then exactly `nsz`. The fixture above
  draws about 3 million times per run, so it took roughly one run in 700; a
  genome-scale MSigDB collection draws enough to hit it far more often than
  that. `Random()` returns a value strictly below 1, and the index is clamped
  besides.

- Count summaries handle singleton and zero-total count vectors with a Gini
  index of zero, avoiding a division by zero or negative QC values (#48).

- Automatic UMI detection now includes the final read position and closes a
  variable region at the read end, allowing terminal UMIs to be counted (#34).

- Three typos in help text users see: `beta_0` was spelled `bata_0`
  (`mle -b/--beta-labels`), `integer` as `interger` (`mle --permutation-round`),
  and `permutation` as `perumtation` (`pathway --permutation`). Also the grammar
  in `mle --day0-label`, `--cnv-est` (both `test` and `mle`) and
  `mle --max-sgrnapergene-permutation`.
- `test -k/--count-table` and `plot -k/--count-table` described the count table as
  an alternative to "sam files". Neither subcommand accepts SAM input; the clause
  was left over from `count`.
- `--control-sgrna`/`--control-gene` claimed, under every subcommand that accepts
  them, to generate "the null distribution of RRA". Each subcommand now describes
  what it actually does with the list. Under `test` the original wording was
  right. Under `mle` the list is used for normalization and to build the
  background gene set the permutation draws on (`mlemultiprocessing.py:245`),
  which is not RRA — though not under `--no-permutation-by-group`, which takes a
  different path (`mlemageck.py:232`) and permutes against all genes. Under
  `plot` it is used for normalization alone. Under `count` it is
  used for normalization, and — when `--day0-label` is given and the `--gmt-file`
  pathway overlaps the library — the quality control step shells out to `test`
  with the same list (`mageckCountQC.py:156`), where it does become RRA's control
  (`crisprFunction.py:561`).
- `test --remove-zero-threshold` said an sgRNA counts as zero when its normalized
  mean falls "below" the threshold. The filter is `x > threshold`, so a mean equal
  to the threshold is removed as well; it now reads "at or below". This matters at
  the default of 0, where it decides whether an exactly-zero sgRNA is dropped.

## [0.3.0] - 2026-08-31

This release changes two command-line behaviors. `--pairguide auto` is no longer
accepted, and `mle` now refuses `-p`/`--PPI-prior`, `-w`/`--PPI-weighting` and
`-e`/`--negative-control` instead of ignoring them. Both previously ran to
completion and produced output that did not reflect the option given; see
**Removed** and **Fixed** below.

### Added

- Bundled the `mageckGSEA` C++ helper (vendored from legacy MAGeCK, BSD
  licensed). It is compiled and installed onto `PATH` alongside `RRA`, so
  `pathway --method gsea` — the default — now works from a pip/source install.
- Pathway smoke tests (`pathway --method gsea` and `--method rra`, end to end)
  and a small GMT test fixture.
- `tests/conftest.py`, which puts this repository's freshly built `RRA` and
  `mageckGSEA` first on `PATH` when they exist. The suite invokes the helpers by
  bare name, so under an editable install — which never runs the compile step —
  it silently tested whichever build came first, including a MAGeCK v1 package's
  identically named binaries. That reported four `mageckGSEA` failures against
  fixes present in the checkout, with nothing in the output naming the binary as
  the cause. The override applies in CI too, since `pip install .` leaves the
  build outputs in the source tree alongside the copies `data_files` installs;
  so the two tests that assert the helpers are *installed* now look at the PATH
  as it stood before the override, keeping a broken installation from passing.
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
- `count` validates the extraction window for `--pairguide` and `--umi` once, up
  front, and refuses to start without it. See issue #29.
- Paired-guide tests for a misplaced `--pg-start-2`/`--pg-end-2` window, for
  `--reverse-complement-2` in both orientations, for the `--pg-min-read`
  threshold including its boundary value, and for `--pg-pair-only` omitted.
  See issue #29.

### Removed

- `--pairguide auto` is no longer accepted, pending a reimplementation. It reused
  the UMI variable-region search, which looks for a run of columns whose base
  composition is near-random. That describes a UMI, not a guide drawn from a
  fixed library, so the window it returned was truncated or never closed, and
  every layout tested produced a header-only `pg_count.txt` with a zero exit
  code. Specify `--pg-start`/`--pg-end` or `--pg-start-2`/`--pg-end-2` instead.
  See issue #29 and issue #32.

### Fixed

- `mle` no longer accepts `-p`/`--PPI-prior`, `-w`/`--PPI-weighting` or
  `-e`/`--negative-control` and then ignores them. All three are read only
  inside the experimental Bayes module, which is disabled, so a run passing them
  exited 0 and wrote a complete result set produced without them — a "PPI on"
  and a "PPI off" run gave the same model, differing only in unseeded
  permutation noise. They are now refused together, exiting non-zero and writing
  nothing, matching what `--bayes` already did. `-e` additionally points at
  `--norm-method control` with `--control-gene`, which is what a user reaching
  for it usually wants — and warns off `--control-sgrna`, which takes sgRNA IDs
  rather than the gene names `-e` accepted. See issue #36.
- `--pairguide firstpair|secondpair` and `--umi firstpair|secondpair` no longer
  accept a missing extraction window. The window was checked inside the per-read
  loop, which logged an error for every read — one line per read on a real FASTQ
  — then sliced an empty string anyway and exited 0 with an empty table. The
  check now runs once, before any FASTQ is opened, and exits non-zero. See
  issue #29.
- `--umi auto` no longer accepts a failed search as a result. A second-read
  search that set only one of the two bounds was treated as success, so the
  failure value `(0, -1)` sliced each read as `[0:-1]` — the whole read but its
  last base — and every read produced a distinct "UMI", filling `umi_count.txt`
  with plausible-looking rows built from unrelated sequence. See issue #29.
- `--umi auto` no longer dies with an uncaught `ValueError` from `max()` on an
  empty sequence when nothing follows the guide on read 1; it reports that there
  is no region to search. See issue #29.
- In `--pairguide secondpair` mode, the second-read window is no longer recorded
  for read pairs whose first read matched no guide. Those were accumulated under
  a `None` key that nothing ever reads back — `mageckcount_printpgdict` skips it
  — so the memory retained grew with the *unmapped* fraction of the FASTQ rather
  than with the library. A run of 100 matched against 5,000 unmatched reads held
  5,000 such entries against 100 useful ones. See issue #29.
- Removed an unreachable `--pairguide` branch in `mageckcount_processonefile`.
  `args.umi` is assigned from `args.pairguide` before it runs, so the preceding
  test always won; had it run, it set `args.umi='none'` and disabled paired-guide
  extraction entirely.
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

[Unreleased]: https://github.com/davidliwei/mageck2/compare/v0.3.0...HEAD
[0.3.0]: https://github.com/davidliwei/mageck2/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/davidliwei/mageck2/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/davidliwei/mageck2/releases/tag/v0.1.0
