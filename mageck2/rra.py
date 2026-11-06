"""Robust Rank Aggregation, in Python.

This is the port of the bundled ``RRA`` C++ helper. It exists so that ranking
can be called in process -- arrays in, an object out -- instead of through
``subprocess`` with a temporary file on each side, which is what
``crisprFunction.py`` does today and what prevents the public API from
returning anything other than a file path.

The port is held to agreement with the C++ rather than to its own internal
consistency: ``tests/parity/`` records goldens from the compiled helper and
``tests/test_rra_parity.py`` compares this module against them, gene by gene.
Several choices here look gratuitously literal and are deliberate, because
agreement is the requirement:

* ``nearest_index`` transcribes ``bTreeSearchingF`` (``math_api.cpp:29``)
  including the guards it re-applies at every recursion level, rather than
  calling ``np.searchsorted``. The C++ function is a nearest-element search,
  not a bisection: it returns the index of the closest element that is
  actually present. ``searchsorted`` returns an insertion point, which differs
  on essentially every value -- measured at 100% of values on the parity
  corpus, continuous data included -- so it is not a near-equivalent that
  happens to disagree on ties. It is a different function.
* ``LehmerStream`` reproduces ``rngs.cpp`` exactly, including the stream the
  helper plants at startup, so the permutation null is not merely drawn from
  the same distribution but is the same sequence of numbers. That turns the
  comparison against the goldens into an exact test rather than a statistical
  one, which is the only kind of test that can be run per gene.
* ``_lo_values`` keeps ``ComputeLoValue``'s quirk of always evaluating the
  smallest percentile even when it already exceeds the alpha cutoff, while
  not counting it as a "good" sgRNA.

The incomplete beta function is the one place the C++ is *not* reproduced
bit for bit, because it does not need to be:
``BetaNoncentralCdf(a, b, 0.0, x, 1e-10)`` is identically
``scipy.special.betainc(a, b, x)``. The noncentrality parameter is 0.0 at every
call site in the C++, so its correction loop never runs. Over a 20,790-point
grid the two agree to 2.1e-11 absolute, and SciPy is the more accurate of the
two: ``betainc(1, 1, x)`` deviates from ``x`` by up to 1.0e-11 in the C++ and
by 0.0 in SciPy.
"""

import os

import re

import numpy as np
from scipy.special import betainc

#: Lehmer generator constants from ``rngs.cpp``. The helper plants
#: ``PlantSeeds(123456)`` on stream 0 at the start of every permutation run
#: (``RRA.cpp:722,877``), which is why its output is reproducible at all.
MODULUS = 2147483647
MULTIPLIER = 48271
DEFAULT_SEED = 123456

#: ``MAX_NAME_LEN`` in ``classdef.h``; names are stored truncated to one less
#: than this, and control-list lookups compare against the truncated form.
MAX_NAME_LEN = 512

#: Defaults mirroring the C++ command line.
DEFAULT_ALPHA = 0.25
DEFAULT_PERMUTATION = 100
DEFAULT_MAX_SGRNA_PER_GENE_PERMUTATION = 100
#: Lower bound on the percentiles that enter a gene's lo-value -- the
#: helper's ``-P``. The default is inert: percentiles lie in [0, 1], so a
#: bound of -1.0 excludes nothing, which is why the option can be absent
#: from every recorded golden without any of them being wrong.
DEFAULT_MIN_PERCENTILE = -1.0

#: Below this many total permutations the helper silently raises the count
#: (``RRA.cpp:281``) unless ``--permutation`` was given explicitly.
AUTO_PERMUTATION_TARGET = 100000

#: Elements of the permutation null to materialize at a time. The null has
#: (permutation + 1) * (number of genes) rows, so a 20,000-gene screen builds
#: 2.02 million of them per sgRNA-count group, and the lo-value computation
#: holds several arrays of the full (rows x sgRNAs) shape at once. Each row is
#: independent, so chunking changes nothing about the result and caps the peak.
#:
#: The bound counts elements rather than rows on purpose. A row limit is not a
#: memory limit: it is only one if every row has the same width, and gene sizes
#: in the wild range from 1 sgRNA to the few hundred that a tiling or saturation
#: library puts on one gene. At 2**18 rows a 100-sgRNA gene would allocate a
#: 262,144 x 100 float64 array -- 210 MB, times the four arrays of that shape
#: the reduction holds -- so the limit would be defeated by exactly the input
#: that needs it. 2**18 elements is about 2 MB per array.
NULL_CHUNK_ELEMENTS = 1 << 18

#: ``strtod`` and ``atoi`` accept a leading numeric prefix and yield 0 when
#: there is none; these recognise that prefix so a malformed optional column
#: is read the way the C++ reads it instead of raising.
_LEADING_FLOAT = re.compile(r"[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?")
_LEADING_INT = re.compile(r"[+-]?\d+")

#: The epsilon the helper brackets each value with when converting a rank to a
#: percentile, and a lo-value to a p value.
EPS = 0.000000001


def name_as_stored(name):
    """Truncate a name the way the C++ stores it, so lookups compare alike."""
    return name[:MAX_NAME_LEN - 1]


def _atof(text):
    """``strtod``'s leading-numeric-prefix rule: no prefix reads as 0.0."""
    match = _LEADING_FLOAT.match(text.strip())
    return float(match.group(0)) if match else 0.0


def _atoi(text):
    """``atoi``'s leading-integer-prefix rule: no prefix reads as 0."""
    match = _LEADING_INT.match(text.strip())
    return int(match.group(0)) if match else 0


class LehmerStream:
    """The ``rngs.cpp`` generator, drawn in vectorized blocks.

    The recurrence is ``state <- MULTIPLIER * state mod MODULUS``; the C++
    evaluates it with Schrage's trick to stay inside a 32-bit signed integer,
    which is an implementation detail of avoiding overflow and does not change
    the sequence.

    Drawing one value at a time in Python would make the permutation null
    slower than the C++ it replaces, so values are produced a block at a time
    by jumping ahead: the k-th state after ``s`` is ``s * MULTIPLIER**k mod
    MODULUS``, and the powers are precomputed once. Both factors are below
    2**31, so each product fits in int64 and the arithmetic is exact -- no
    floating point enters until the final division.
    """

    def __init__(self, seed=DEFAULT_SEED, block=1 << 16):
        if not 1 <= seed <= MODULUS - 1:
            raise ValueError(
                "seed must be between 1 and %d; got %r" % (MODULUS - 1, seed))
        self.state = int(seed)
        self._block = int(block)
        powers = np.empty(self._block, dtype=np.int64)
        current = 1
        for i in range(self._block):
            current = (current * MULTIPLIER) % MODULUS
            powers[i] = current
        self._powers = powers

    def uniform(self, n):
        """Return the next ``n`` values, each in (0, 1)."""
        out = np.empty(n, dtype=np.float64)
        taken = 0
        while taken < n:
            size = min(self._block, n - taken)
            states = (self.state * self._powers[:size]) % MODULUS
            out[taken:taken + size] = states
            self.state = int(states[size - 1])
            taken += size
        return out / MODULUS


def nearest_index(sorted_values, queries):
    """Index of the element of ``sorted_values`` nearest to each query.

    A vectorized transcription of ``bTreeSearchingF``. The loop below performs
    the same three tests in the same order as each recursive call -- clamp low,
    clamp high, then resolve a two-element window by distance with ties going
    to the upper index -- before bisecting. Running them in that order matters:
    the clamps are what make an exact hit return the index it was found at,
    which is the function's tie behaviour.
    """
    a = np.asarray(sorted_values, dtype=np.float64)
    v = np.atleast_1d(np.asarray(queries, dtype=np.float64))
    n = a.shape[0]
    if n == 0:
        raise ValueError("cannot search an empty array")
    if n == 1:
        return np.zeros(v.shape, dtype=np.int64)

    lo = np.zeros(v.shape, dtype=np.int64)
    hi = np.full(v.shape, n - 1, dtype=np.int64)
    out = np.empty(v.shape, dtype=np.int64)
    pending = np.ones(v.shape, dtype=bool)

    while True:
        hit = pending & (v <= a[lo])
        out[hit] = lo[hit]
        pending &= ~hit

        hit = pending & (v >= a[hi])
        out[hit] = hi[hit]
        pending &= ~hit

        hit = pending & (hi - lo <= 1)
        if hit.any():
            take_lo = np.abs(a[hi[hit]] - v[hit]) > np.abs(a[lo[hit]] - v[hit])
            out[hit] = np.where(take_lo, lo[hit], hi[hit])
            pending &= ~hit

        if not pending.any():
            return out

        mid = (lo + hi) // 2
        upper = pending & (v >= a[mid])
        lower = pending & ~upper
        lo = np.where(upper, mid, lo)
        hi = np.where(lower, mid, hi)


def percentiles(sorted_values, values):
    """Percentile of each value within its sorted list.

    ``ProcessGroups`` (``RRA.cpp:414-418``) brackets each value by EPS and
    averages the two indices, which spreads a run of tied values across the
    midpoint of the block they occupy instead of putting them all at one end.
    """
    n = len(sorted_values)
    v = np.asarray(values, dtype=np.float64)
    low = nearest_index(sorted_values, v - EPS)
    high = nearest_index(sorted_values, v + EPS)
    return (low + high + 1.0) / (n * 2.0)


def _lo_values(matrix, alpha, min_percentile=DEFAULT_MIN_PERCENTILE):
    """Lo-values and good-sgRNA counts for each row of percentiles.

    ``ComputeLoValue`` (``RRA.cpp:470``) walks a gene's percentiles in
    ascending order and stops at the first one outside ``[min_percentile,
    alpha]`` -- except at index 0, which is evaluated regardless so that a gene
    whose best sgRNA already misses the cutoff still receives a lo-value rather
    than the sentinel 1.0. An index evaluated that way does not count towards
    ``goodsgrna``, so a gene can have a lo-value below 1 and no good sgRNAs at
    all.

    Because the stop is at the *first* out-of-range index, the evaluated set is
    always a prefix, which is what makes this expressible without a loop.
    """
    x = np.sort(np.asarray(matrix, dtype=np.float64), axis=1)
    rows, k = x.shape
    if k == 0:
        return np.ones(rows), np.zeros(rows, dtype=np.int64)

    in_range = (x <= alpha) & (x >= min_percentile)
    # The prefix ends at the first out-of-range index at or after 1; index 0 is
    # evaluated either way.
    if k > 1:
        out_of_range = ~in_range[:, 1:]
        any_bad = out_of_range.any(axis=1)
        first_bad = np.where(any_bad, out_of_range.argmax(axis=1) + 1, k)
    else:
        first_bad = np.full(rows, k)
    evaluated = np.arange(k)[None, :] < first_bad[:, None]

    a = np.arange(1, k + 1, dtype=np.float64)
    b = (k - np.arange(k)).astype(np.float64)
    cdf = betainc(a[None, :], b[None, :], x)

    lo = np.where(evaluated, cdf, np.inf).min(axis=1)
    lo = np.minimum(lo, 1.0)
    good = (in_range & evaluated).sum(axis=1).astype(np.int64)
    return lo, good


def _benjamini_hochberg(p):
    """``calculateFDRVector`` (``RRA.cpp:811``).

    Standard Benjamini-Hochberg with one wrinkle worth keeping: only the
    largest adjusted value is clamped to 1.0, and the rest are brought below it
    by the monotonicity pass rather than by being clamped individually.
    """
    p = np.asarray(p, dtype=np.float64)
    n = p.shape[0]
    if n == 0:
        return np.zeros(0)
    order = np.argsort(p, kind="stable")
    adjusted = p[order] * n / np.arange(1, n + 1, dtype=np.float64)
    adjusted[n - 1] = min(adjusted[n - 1], 1.0)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    out = np.empty(n, dtype=np.float64)
    out[order] = adjusted
    return out


class RRAResult:
    """One row per gene, in the order the C++ writes them.

    ``to_text`` renders the historical file format, which is what lets goldens
    recorded from the C++ serve as the oracle for this implementation.
    """

    HEADER = "group_id\titems_in_group\tlo_value\tp\tFDR\tgoodsgrna"

    def __init__(self, names, items, lo_value, p, fdr, goodsgrna, excluded,
                 suppressed, missing_control=()):
        self.names = list(names)
        self.items = np.asarray(items, dtype=np.int64)
        self.lo_value = np.asarray(lo_value, dtype=np.float64)
        self.p = np.asarray(p, dtype=np.float64)
        self.fdr = np.asarray(fdr, dtype=np.float64)
        self.goodsgrna = np.asarray(goodsgrna, dtype=np.int64)
        #: Genes held out of FDR adjustment by the --min-*-goodsgrna filters.
        #: Their FDR prints as NA rather than as a number.
        self.excluded = np.asarray(excluded, dtype=bool)
        #: Genes made up entirely of control sgRNAs, which the C++ omits from
        #: the output altogether (``fileio.cpp:331``).
        self.suppressed = np.asarray(suppressed, dtype=bool)
        #: Control sgRNA names that never appeared in the ranked list, in the
        #: order the C++ reports them. They contribute nothing to the null,
        #: so this changes no number -- but a user who mistypes a name, or
        #: whose control guides were filtered out upstream, is told. The C++
        #: warns here (``RRA.cpp:460``) and callers surface it; dropping the
        #: list silently would make the Python backend quieter than the
        #: helper it replaces about something worth knowing.
        self.missing_control = list(missing_control)

    def __len__(self):
        return len(self.names)

    def to_text(self):
        lines = [self.HEADER]
        order = np.argsort(self.lo_value, kind="stable")
        for i in order:
            if self.suppressed[i]:
                continue
            fdr = "NA" if self.excluded[i] else "%f" % self.fdr[i]
            lines.append("%s\t%d\t%10.4e\t%10.4e\t%s\t%d" % (
                self.names[i], self.items[i], self.lo_value[i], self.p[i],
                fdr, self.goodsgrna[i]))
        return "\n".join(lines) + "\n"


def _resolve_permutation_count(requested, ngenes, max_per_group,
                               user_specified):
    """Reproduce the two adjustments the C++ makes to the permutation count.

    ``RRA.cpp:281`` quietly raises the count when the screen is small enough
    that the default would not buy usable p value resolution, and
    ``RRA.cpp:293`` lowers it when ``--max-permutation-per-group`` asks for a
    bound. Both are skipped when ``--permutation`` was given explicitly, so a
    caller who states a count gets that count.

    The kernels run ``(requested + 1) * ngenes`` permutations per group, not
    ``requested * ngenes``, which is why the bound is compared against the
    former.
    """
    count = requested
    if not user_specified and count * ngenes < AUTO_PERMUTATION_TARGET:
        count = AUTO_PERMUTATION_TARGET // ngenes + 1
    if (max_per_group > 0 and not user_specified
            and (count + 1) * ngenes > max_per_group):
        count = max_per_group // ngenes - 1
        if count < 1:
            count = 1
    return count


def _control_percentile_array(control_names, percentile_by_sgrna):
    """The array the permutation resamples control percentiles from.

    Order matters and is not arbitrary: the C++ walks ``ControlSeqMap``, a
    ``std::map<string,int>``, so the array is built in lexicographic order of
    sgRNA name, and ``rand_ctl_index`` indexes into that order. Building it in
    any other order would draw from the same multiset but in a different
    sequence, and the comparison against the goldens would stop being exact.

    Control sgRNAs absent from the ranked list are skipped, which is what the
    C++ warns about and then ignores.
    """
    values = []
    missing = []
    for name in sorted(set(control_names)):
        if name in percentile_by_sgrna:
            values.append(percentile_by_sgrna[name])
        else:
            missing.append(name)
    return np.asarray(values, dtype=np.float64), missing


def rank(items, groups, values, lists=None, alpha=DEFAULT_ALPHA, control=None,
         permutation=DEFAULT_PERMUTATION, permutation_specified=False,
         by_group=True, skip_genes=(), min_number_goodsgrna=0,
         min_percentage_goodsgrna=0.0,
         max_sgrna_per_gene_permutation=DEFAULT_MAX_SGRNA_PER_GENE_PERMUTATION,
         max_permutation_per_group=0, seed=DEFAULT_SEED,
         # Last and keyword-only, although it belongs next to ``alpha``: the
         # two are the helper's -p and -P. Putting it there rebinds every
         # positional argument after it, so ``rank(items, groups, values,
         # lists, alpha, control)`` would pass a control list in as a
         # percentile bound. Keyword-only also means the next option added
         # here cannot repeat that.
         *, min_percentile=DEFAULT_MIN_PERCENTILE):
    """Rank groups by robust rank aggregation.

    ``items`` are sgRNA names, ``groups`` the gene each belongs to, ``values``
    the statistic to rank on (smaller ranks first), and ``lists`` the pool each
    sgRNA was measured in -- percentiles are computed within a pool. All four
    are parallel sequences.

    ``alpha`` and ``min_percentile`` are the helper's ``-p`` and ``-P``: a
    gene's lo-value is computed over the prefix of its percentiles that lies
    within ``[min_percentile, alpha]``. The lower bound applies to the
    permutation null as well as to the observed lo-values, because the C++
    passes ``minPercentile`` into ``ComputePermutation`` and
    ``ComputePermutationbyNumSG`` too (``RRA.cpp:306-308``). Bounding only the
    observed side would test real genes against a null built under a
    different rule.

    Returns an :class:`RRAResult`.
    """
    items = [name_as_stored(str(x)) for x in items]
    # Group fields are kept whole: each is split on "," below and every
    # sub-name is a grouping key in full. Truncation applies to the stored
    # label only, never to the key.
    groups = [str(x) for x in groups]
    values = np.asarray(values, dtype=np.float64)
    if lists is None:
        lists = np.zeros(len(items), dtype=np.int64)
    lists = np.asarray(lists)
    if not (len(items) == len(groups) == len(values) == len(lists)):
        raise ValueError("items, groups, values and lists must be parallel")
    if len(items) == 0:
        raise ValueError("no sgRNAs to rank")
    if min_percentile > 1.0:
        # The helper refuses this rather than running (RRA.cpp:246). Every
        # percentile is at most 1, so a bound above 1 excludes all of them and
        # would report a lo-value of 1.0 for every gene -- an answer that looks
        # like a result and is not one.
        raise ValueError(
            "minimum percentile (-P) must be at most 1.0; got %r"
            % (min_percentile,))

    # Percentiles are assigned within each pool, against that pool's sorted
    # values -- not globally.
    pct = np.empty(len(values), dtype=np.float64)
    for pool in np.unique(lists):
        mask = lists == pool
        pct[mask] = percentiles(np.sort(values[mask]), values[mask])

    # A row's group field is a comma-separated list of every group the row
    # belongs to, and the row is a member of all of them -- the C++ splits on
    # "," and loops over the parts, both when counting groups and when filling
    # them (fileio.cpp:107-141 and 266-287). This is not a corner case:
    # mageck2's pathway mode depends on it, because pathwayFunc.py:155-162
    # joins every pathway containing a gene into one field, so a gene in three
    # pathways arrives as "P1,P2,P3" and has to count in all three. Treating
    # the field as one name would invent a pathway that does not exist and
    # leave the three real ones a member short each.
    #
    # A repeated sub-name is not collapsed: the C++ increments itemNum once per
    # occurrence and appends the row once per occurrence, so "P1,P1" really
    # does give P1 the row twice.
    #
    # Keys are the full sub-names. Only the label kept for output is
    # truncated, which is what copyName does (fileio.cpp:44-58) and what its
    # warning promises -- "Grouping is unaffected -- it keys on the full name".
    # Two sub-names agreeing to MAX_NAME_LEN-1 therefore stay two groups that
    # happen to print the same label, rather than merging into one.
    #
    # Groups keep the order of first appearance, scanning rows in order and
    # each row's sub-names left to right: that is the order the C++ assigns
    # tmpGroupNum in, and therefore the order its pooled permutation consumes
    # random numbers in.
    gene_index = {}
    gene_names = []
    row_groups = []
    for g in groups:
        indices = []
        for sub in g.split(","):
            if sub not in gene_index:
                gene_index[sub] = len(gene_names)
                gene_names.append(name_as_stored(sub))
            indices.append(gene_index[sub])
        row_groups.append(indices)
    ngenes = len(gene_names)
    members = [[] for _ in range(ngenes)]
    for i, indices in enumerate(row_groups):
        for gi in indices:
            members[gi].append(i)

    percentile_by_sgrna = {}
    for i, name in enumerate(items):
        percentile_by_sgrna[name] = pct[i]

    control_values = None
    missing_control = []
    control_set = set()
    if control is not None:
        control_set = {name_as_stored(str(c)) for c in control if str(c).strip()}
        control_values, missing_control = _control_percentile_array(
            control_set, percentile_by_sgrna)
        if control_values.size == 0:
            # The C++ computes rand_ctl_index = -1 here and reads
            # control_prob_array[-1] (RRA.cpp:761-762, issue #51). Failing is
            # the deliberate choice: falling back to a uniform null would
            # report p values that are not the control-based ones that were
            # asked for, and say nothing about it.
            raise ValueError(
                "no control sgRNA in the control list was found in the ranked "
                "list, so there is no control distribution to permute against")

    counts = np.array([len(m) for m in members], dtype=np.int64)
    control_counts = np.zeros(ngenes, dtype=np.int64)
    if control_set:
        for gi, rows in enumerate(members):
            control_counts[gi] = sum(1 for r in rows if items[r] in control_set)

    lo_value = np.ones(ngenes, dtype=np.float64)
    goodsgrna = np.zeros(ngenes, dtype=np.int64)
    for k in np.unique(counts):
        gi = np.flatnonzero(counts == k)
        block = np.stack([pct[members[g]] for g in gi])
        lo, good = _lo_values(block, alpha, min_percentile)
        lo_value[gi] = lo
        goodsgrna[gi] = good
        if k <= 1:
            # ComputeLoValue is used rather than ComputeLoValue_Prob whenever
            # a gene has at most one sgRNA; with prob fixed at 1 throughout
            # MAGeCK2 the two agree anyway, so nothing branches here.
            pass

    skip = {name_as_stored(str(s)) for s in skip_genes}
    # A gene made up entirely of control sgRNAs is skipped and suppressed.
    suppressed = control_counts >= counts
    for gi in np.flatnonzero(suppressed):
        skip.add(gene_names[gi])

    n_pass = _resolve_permutation_count(
        permutation, ngenes, max_permutation_per_group, permutation_specified)
    scan_pass = n_pass + 1

    rng = LehmerStream(seed)
    p = np.ones(ngenes, dtype=np.float64)
    skipped = np.array([g in skip for g in gene_names], dtype=bool)

    if by_group:
        strata = sorted({int(k) for k, s in zip(counts, skipped) if not s})
        null = None
        for k in strata:
            if k <= max_sgrna_per_gene_permutation:
                null = np.sort(_null_lo_values(
                    rng, scan_pass * ngenes, k, alpha, control_values,
                    min_percentile=min_percentile))
            if null is None:
                continue
            target = (counts == k) & ~skipped
            if target.any():
                p[target] = _p_from_null(null, lo_value[target])
    else:
        # The pooled kernel draws for every gene in input order, including
        # skipped ones, so the null mixes the library's actual sgRNA counts.
        per_pass = int(counts.sum())
        offsets = np.concatenate([[0], np.cumsum(counts)])
        lo_null = np.empty(scan_pass * ngenes, dtype=np.float64)
        # Chunked over passes rather than rows, because one pass draws for
        # every gene in turn and cannot be split without disturbing the order
        # the generator is consumed in.
        # per_pass is the width of one pass over every gene, so dividing by
        # it already bounds elements rather than rows.
        per_chunk = max(1, NULL_CHUNK_ELEMENTS // max(per_pass, 1))
        first = 0
        while first < scan_pass:
            npass = min(per_chunk, scan_pass - first)
            draws = rng.uniform(npass * per_pass).reshape(npass, per_pass)
            for k in np.unique(counts):
                gi = np.flatnonzero(counts == k)
                block = np.concatenate(
                    [draws[:, offsets[g]:offsets[g] + k] for g in gi], axis=0)
                if control_values is not None:
                    block = _resample_control(block, control_values)
                lo, _ = _lo_values(block, alpha, min_percentile)
                # block stacks each gene's npass rows in turn, so reshaping to
                # (gene, pass) recovers which lo-value belongs where.
                where = (first + np.arange(npass))[None, :] * ngenes + gi[:, None]
                lo_null[where] = lo.reshape(len(gi), npass)
            first += npass
        null = np.sort(lo_null)
        p = _p_from_null(null, lo_value)

    # adjustFDR (RRA.cpp:641): genes with too few good sgRNAs are held out of
    # the adjustment entirely rather than adjusted and then ignored.
    with np.errstate(invalid="ignore", divide="ignore"):
        fraction = np.where(counts > 0, goodsgrna / counts, 0.0)
    excluded = (fraction < min_percentage_goodsgrna) | (
        goodsgrna < min_number_goodsgrna)
    fdr = np.ones(ngenes, dtype=np.float64)
    keep = ~excluded
    if keep.any():
        fdr[keep] = _benjamini_hochberg(p[keep])

    return RRAResult(gene_names, counts, lo_value, p, fdr, goodsgrna,
                     excluded, suppressed, missing_control=missing_control)


def _null_lo_values(rng, rows, k, alpha, control_values,
                    min_percentile=DEFAULT_MIN_PERCENTILE,
                    chunk=NULL_CHUNK_ELEMENTS):
    """Lo-values for ``rows`` permuted genes of ``k`` sgRNAs each.

    Drawn and reduced in chunks of at most ``chunk`` elements -- so at most
    ``chunk // k`` rows, never ``chunk`` rows of unbounded width. The generator
    is consumed strictly in order, so chunking does not disturb the sequence,
    and each row's lo-value depends only on its own draws: the chunk boundary
    is invisible in the result.
    """
    out = np.empty(rows, dtype=np.float64)
    per_chunk = max(1, chunk // max(k, 1))
    done = 0
    while done < rows:
        take = min(per_chunk, rows - done)
        block = rng.uniform(take * k).reshape(take, k)
        if control_values is not None:
            block = _resample_control(block, control_values)
        out[done:done + take] = _lo_values(block, alpha, min_percentile)[0]
        done += take
    return out


def _resample_control(block, control_values):
    """Replace uniform draws with resampled control-sgRNA percentiles.

    The clamp mirrors ``RRA.cpp:762``. It cannot fire here -- ``Random()``
    returns strictly below 1 because its state stays in ``1..MODULUS-1`` -- and
    is kept because the index arithmetic is the C++'s, not because the bound is
    reachable.
    """
    n = control_values.shape[0]
    idx = (block * n).astype(np.int64)
    np.minimum(idx, n - 1, out=idx)
    return control_values[idx]


def _p_from_null(null, lo_values):
    """Empirical p value of each lo-value against the sorted null."""
    low = nearest_index(null, lo_values - EPS)
    high = nearest_index(null, lo_values + EPS)
    return (low + high + 1.0) / 2.0 / null.shape[0]


#: Flags accepted by ``rank_from_file``, mapping the C++ option spelling to
#: the keyword it sets. Kept as data so the set of supported options is
#: visible in one place rather than spread through a parser.
_VALUE_FLAGS = {
    "-P": ("min_percentile", float),
    "--permutation": ("permutation", int),
    "--min-number-goodsgrna": ("min_number_goodsgrna", int),
    "--min-percentage-goodsgrna": ("min_percentage_goodsgrna", float),
    "--max-sgrnapergene-permutation": ("max_sgrna_per_gene_permutation", int),
    "--max-permutation-per-group": ("max_permutation_per_group", int),
}


def _parse_flags(flags):
    """Translate a list of C++ command-line flags into ``rank`` keywords."""
    kwargs = {}
    skip = []
    i = 0
    flags = list(flags or ())
    while i < len(flags):
        flag = flags[i]
        if flag == "--no-permutation-by-group":
            kwargs["by_group"] = False
            i += 1
        elif flag == "--skip-gene":
            skip.append(flags[i + 1])
            i += 2
        elif flag in _VALUE_FLAGS:
            key, cast = _VALUE_FLAGS[flag]
            kwargs[key] = cast(flags[i + 1])
            if flag == "--permutation":
                # The auto-raise and the ceiling are both suppressed when the
                # count was asked for explicitly.
                kwargs["permutation_specified"] = True
            i += 2
        else:
            raise ValueError("unsupported RRA flag: %r" % (flag,))
    if skip:
        kwargs["skip_genes"] = skip
    return kwargs


def read_input(path):
    """Read the helper's input format.

    ``<item> <group> <list> <value> [<prob>] [<iscounted>]``, behind a header
    row that is discarded -- which is what ``ReadFile`` does.

    The two optional columns are parsed and checked rather than skipped over.
    mageck2 writes them itself (``crisprFunction.py:505`` and ``:529``) but
    always with the values that make them inert: ``prob`` is the literal
    ``1``, and rows are written only when ``chosen`` is 1, so no file this
    package produces is refused here. A file that does carry a meaningful
    value is a file whose statistics this port cannot reproduce, and refusing
    it is the point -- read silently, those columns change the answer:

    * ``prob`` below 1 sends the C++ into ``ComputeLoValue_Prob``
      (``RRA.cpp:525``), which sums over all ``2**n`` subsets of a gene's
      sgRNAs. Besides being unimplemented here, it is unusable above a few
      dozen sgRNAs per gene, and it pairs the probabilities with the
      *unsorted* percentiles while reducing over the sorted ones, so a port
      would have to reproduce that pairing to agree with it.
    * ``chosen`` of 0 drops the row from its pool's percentile distribution
      and from its gene's lo-value, while still counting it in the gene's
      sgRNA total -- so the lo-value is computed over fewer sgRNAs than the
      total it is reported against, and the permutation null for the gene is
      built from the smaller count (``RRA.cpp:412-413``, ``:646``,
      ``:908-918``).

    Supporting either means recording goldens for it while the C++ is still
    here to record them from, which is Phase 3b work at the earliest.
    """
    items, groups, pools, values = [], [], [], []
    with open(path) as handle:
        for lineno, line in enumerate(handle, 1):
            if lineno == 1 or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                fields = line.split()
            if len(fields) < 4:
                raise ValueError(
                    "%s:%d: expected 4 columns, got %d"
                    % (path, lineno, len(fields)))
            items.append(fields[0])
            groups.append(fields[1])
            pools.append(fields[2])
            values.append(float(fields[3]))
            # atof/atoi, like the C++: a field that is not a number reads as
            # zero rather than raising, so a non-numeric prob is caught by the
            # value check below rather than by float().
            if len(fields) > 4 and _atof(fields[4]) != 1.0:
                raise ValueError(
                    "%s:%d: sgRNA probabilities other than 1 are not "
                    "supported by the Python RRA backend (got %r). Run with "
                    "MAGECK2_RRA_IMPL=cpp, or drop the column."
                    % (path, lineno, fields[4]))
            if len(fields) > 5 and _atoi(fields[5]) != 1:
                raise ValueError(
                    "%s:%d: unchosen sgRNAs (column 6 other than 1) are not "
                    "supported by the Python RRA backend (got %r). Run with "
                    "MAGECK2_RRA_IMPL=cpp, or remove the row."
                    % (path, lineno, fields[5]))
    if not items:
        raise ValueError("%s contains no data rows" % path)
    pool_ids = {name: i for i, name in enumerate(dict.fromkeys(pools))}
    return (items, groups, np.asarray(values, dtype=np.float64),
            np.asarray([pool_ids[p] for p in pools], dtype=np.int64))


def read_control(path):
    """Read a control sgRNA list, one name per line.

    Blank lines are dropped. The C++ reader does not drop them: it loops on
    ``while (!fh.eof())`` and so inserts an empty-named entry for the newline
    terminating the last line, which inflates its reported count by one and
    makes it warn about an sgRNA with no name. Dropping them here is a
    deliberate divergence in a diagnostic message, not in any p value -- an
    empty name matches no sgRNA and so never entered the null.
    """
    with open(path) as handle:
        return [line.strip() for line in handle if line.strip()]


def rank_from_file(path, alpha=DEFAULT_ALPHA, control=None, flags=None,
                   seed=DEFAULT_SEED):
    """Run RRA over an input file, mirroring the compiled helper's CLI.

    ``control`` is a list of sgRNA names, or a path to a file of them.
    """
    items, groups, values, pools = read_input(path)
    if isinstance(control, str):
        control = read_control(control)
    kwargs = _parse_flags(flags)
    return rank(items, groups, values, lists=pools, alpha=alpha,
                control=control, seed=seed, **kwargs)
