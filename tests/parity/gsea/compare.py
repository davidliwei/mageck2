"""Column-aware comparison of two ``mageckGSEA`` outputs.

Comparing these files byte for byte would be wrong in both directions. It is
too strict, because every number is written through ``ostream <<`` with default
flags -- six significant digits -- so the file already discards precision the
implementations may legitimately differ in below. And it is too weak, because
a byte diff says nothing about *which* pathways moved, while reporting a
catastrophic difference for a case where two pathways that tie on both sort
keys came out in the other order.

That last point is not hypothetical and is the reason this module exists in the
shape it does. ``testallgseas`` sorts with ``std::sort`` on
``(p_permutation, -ES)``, which is unstable, and pathways do tie on both keys
at once -- three such pairs in the ``tie_extreme`` golden, three in
``big_sets``, two in ``many_pathways``. Measured on a 60-set input at ``-c 0``
while designing this harness: two pathways with bit-identical ES and identical
p came out in opposite orders between two implementations, which shifted the
``Ranking`` column of both and changed nothing else. So:

* ``Ranking`` is never compared directly. What is compared is whether the two
  implementations order each *distinguishable* pair of pathways the same way.
  A pair that ties on both keys carries no ordering claim, and breaking it
  differently is not a disagreement. Nor is the column always an ordering at
  all: without ``-s`` every row gets a literal 0 and the single-pathway branch
  gets a literal 1, so a column of identical values is recognised and skipped
  rather than read as "everything ties for first".
* ``FDR`` *is* compared exactly, because it is not affected. The step-down pass
  takes a running minimum from the bottom, which collapses every tied run to a
  single value regardless of the order within the run.

Three levels of agreement are reported:

1. **Structural** -- the same pathways, with the same ``Size`` and the same
   ``Hits`` (the index at which the enrichment score peaked).
2. **Numeric** -- each column within a tolerance matched to the precision the
   file carries, with the worst observed deviation reported even when it
   passes, so drift that stays inside tolerance is still visible.
3. **Scientific** -- the p-value ranking and the called pathway sets agree.
   This is the level a user would notice.
"""

import dataclasses

import numpy as np

from ..stats import inversions, spearman

#: Columns of the output, in the order ``calculategsea.cpp:488`` writes them.
COLUMNS = ("Pathway", "Size", "ES", "p", "p_permutation", "FDR", "Ranking",
           "Hits", "LFC")

#: ``ostream <<`` with default flags carries six significant digits. Two
#: implementations whose true values round to the same printed digits parse
#: back to *identical* floats, so any non-zero deviation here means the printed
#: digits differ by at least one unit in the sixth significant digit.
#:
#: That unit is worth 1e-5/m relatively, where m is the mantissa -- so between
#: 1e-6 and 1e-5 depending on the value. 2e-5 therefore accepts two such units
#: at the least favourable mantissa and nothing more. (The RRA harness's 1e-4
#: is the same statement one digit coarser.)
RTOL_PRINTED = 2e-5

#: Absolute floor for ``ES``, below which two values are the same number.
#:
#: This one is not about print precision, and it is the most important constant
#: in the file. ``getgseascore_fast`` computes each step as
#: ``decvalue*(ni-i) + sumincval``, which compilers are free to contract into a
#: fused multiply-add -- and clang at ``-O3`` does, on any target where FMA is
#: in the baseline instruction set. An FMA rounds once instead of twice, so it
#: gives a different answer in the last bit.
#:
#: That is normally invisible, because ``gs_ret::gscore`` truncates to a 32-bit
#: float and 1e-17 does not survive the truncation. It becomes visible exactly
#: where the true score is itself near zero: a gene set whose score cancels to
#: 0 is reported as 0 by a non-contracted build and as 5.55112e-17 by a
#: contracted one, and float32 has plenty of resolution to keep the difference.
#:
#: This was measured, not assumed. Rebuilding the same source with
#: ``-ffp-contract=off`` changes three of this corpus's twenty goldens, by
#: 2.776e-17 and 5.551e-17 -- so the helper does not even agree with *itself*
#: across compiler flags, let alone across architectures (FMA is baseline on
#: arm64 and not on baseline x86-64). Requiring a port to match the recorded
#: digits would therefore require it to match one machine's compiler, and
#: would fail the C++ against its own goldens on CI.
#:
#: 1e-15 is an order of magnitude above the largest difference this can
#: produce -- a single contraction error is bounded by half an ulp of a product
#: no larger than 1, so about 1.1e-16 -- and eight orders of magnitude below
#: the float32 resolution of any ES that is not numerically zero. Nothing a
#: port could get genuinely wrong hides under it.
ATOL_ES = 1e-15

#: ``p`` needs no such floor. It is ``BetaNoncentralCdf`` of an integer and an
#: exact ratio of integers, so it depends on the enrichment score only through
#: which index the score peaked at -- and that index was identical across both
#: builds and every case here.

#: ``p_permutation`` is a count of exceedances over a fixed number of
#: permutations, driven by a seeded generator whose consumption order is fully
#: determined. There is no reason for two correct implementations to differ by
#: even one permutation, so it is compared exactly: both tolerances are zero,
#: not just the absolute one. A relative tolerance here would quietly accept a
#: one-permutation difference at large p for a large enough ``-p``.
ATOL_PERMUTATION = 0.0
RTOL_PERMUTATION = 0.0

#: FDR cutoffs at which the called pathway sets are compared.
HIT_THRESHOLDS = (0.05, 0.1, 0.25)

#: Columns compared as exact integers.
INTEGER_COLUMNS = ("Size", "Hits")

#: Columns compared numerically: ``(name, uses_printed_rtol, atol)``. The
#: columns flagged ``True`` are compared at whatever relative tolerance the
#: caller passes; ``p_permutation`` is not, because it is exact.
NUMERIC_COLUMNS = (("ES", True, ATOL_ES), ("p", True, 0.0),
                   ("p_permutation", False, ATOL_PERMUTATION),
                   ("FDR", True, 0.0))


def parse(text):
    """Parse ``mageckGSEA`` output into ``({pathway: row}, row_order)``.

    The pathway name is allowed to be the empty string: ``parsegmt`` stores a
    pathway under that name for a blank GMT line, and the
    ``gmt_blank_line`` corpus case depends on being able to read that row back.
    """
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines:
        raise ValueError("empty mageckGSEA output")
    header = [f.strip() for f in lines[0].split("\t")]
    if tuple(header) != COLUMNS:
        raise ValueError("unexpected mageckGSEA output header: %r" % (header,))

    rows = {}
    order = []
    for ln in lines[1:]:
        f = ln.split("\t")
        if len(f) != len(COLUMNS):
            raise ValueError("malformed mageckGSEA output row: %r" % (ln,))
        name = f[0].strip()
        rows[name] = {
            "Size": int(f[1]),
            "ES": float(f[2]),
            "p": float(f[3]),
            "p_permutation": float(f[4]),
            "FDR": float(f[5]),
            "Ranking": int(f[6]),
            "Hits": int(f[7]),
            "LFC": float(f[8]),
        }
        order.append(name)
    return rows, order


@dataclasses.dataclass
class Deviation:
    """Worst observed disagreement in one numeric column."""

    column: str
    max_abs: float = 0.0
    max_rel: float = 0.0
    worst: str = ""
    n_outside: int = 0

    def __str__(self):
        s = "%-14s max_abs=%.3e max_rel=%.3e" % (
            self.column, self.max_abs, self.max_rel)
        if self.worst:
            s += "  worst=%s" % self.worst
        if self.n_outside:
            s += "  OUTSIDE TOLERANCE: %d pathway(s)" % self.n_outside
        return s


@dataclasses.dataclass
class Report:
    """The outcome of comparing a candidate against a golden."""

    case: str
    missing: list = dataclasses.field(default_factory=list)
    extra: list = dataclasses.field(default_factory=list)
    structural: list = dataclasses.field(default_factory=list)
    deviations: dict = dataclasses.field(default_factory=dict)
    spearman: float = float("nan")
    rank_inversions: int = 0
    ranking_disagreements: int = 0
    tied_pairs: int = 0
    #: False when the output mode writes no ordering into ``Ranking`` -- see
    #: :func:`is_ranked`. The column is then neither checked nor held against
    #: the candidate.
    ranking_checked: bool = True
    order_matches: bool = True
    byte_identical: bool = False
    hits: dict = dataclasses.field(default_factory=dict)

    @property
    def ok(self):
        """True when every asserted level agrees.

        Row order and the ``Ranking`` column are excluded from the numeric
        checks and enter only through ``ranking_disagreements``, which counts
        pairs the golden actually distinguishes. See the module docstring.
        """
        return (not self.missing and not self.extra and not self.structural
                and all(d.n_outside == 0 for d in self.deviations.values())
                and self.ranking_disagreements == 0
                and all(h["only_golden"] == [] and h["only_candidate"] == []
                        for h in self.hits.values()))

    def render(self):
        out = ["case %s: %s" % (self.case, "AGREE" if self.ok else "DISAGREE")]
        if self.missing:
            out.append("  pathways missing from candidate: %d %s"
                       % (len(self.missing), self.missing[:8]))
        if self.extra:
            out.append("  pathways only in candidate: %d %s"
                       % (len(self.extra), self.extra[:8]))
        for line in self.structural[:12]:
            out.append("  " + line)
        if len(self.structural) > 12:
            out.append("  ... %d more structural mismatches"
                       % (len(self.structural) - 12))
        for dev in self.deviations.values():
            out.append("  " + str(dev))
        out.append("  p_permutation ranking: spearman=%.6f inversions=%d"
                   % (self.spearman, self.rank_inversions))
        if self.ranking_checked:
            out.append("  Ranking column: disagreements=%d over "
                       "distinguishable pairs (%d pair(s) tie on both sort "
                       "keys and are exempt)"
                       % (self.ranking_disagreements, self.tied_pairs))
        else:
            out.append("  Ranking column: not an ordering in this output "
                       "mode, so not checked")
        out.append("  row order: %s   byte-identical: %s"
                   % ("same" if self.order_matches else "differs",
                      self.byte_identical))
        for thr, h in sorted(self.hits.items()):
            out.append("  hits at FDR<=%.2f: golden=%d candidate=%d "
                       "only_golden=%s only_candidate=%s"
                       % (thr, h["n_golden"], h["n_candidate"],
                          h["only_golden"][:6], h["only_candidate"][:6]))
        return "\n".join(out)


def _numeric_deviation(column, golden, candidate, names, rtol, atol):
    dev = Deviation(column=column)
    for name in names:
        g = golden[name][column]
        c = candidate[name][column]
        abs_err = abs(g - c)
        if g != 0.0:
            rel_err = abs_err / abs(g)
        else:
            rel_err = 0.0 if c == 0.0 else float("inf")
        if abs_err > dev.max_abs:
            dev.max_abs = abs_err
            dev.worst = name
        dev.max_rel = max(dev.max_rel, rel_err)
        if abs_err > atol + rtol * abs(g):
            dev.n_outside += 1
    return dev


def is_ranked(rows):
    """Whether this output's ``Ranking`` column carries an ordering at all.

    It does not always. ``testallgseas`` writes ``Ranking`` as a literal 0 for
    every row when ``-s`` is absent (``calculategsea.cpp:494``), because
    without the sort there is no rank to report; the single-pathway branch in
    ``main`` writes a literal 1. In both cases every row carries the same
    value, which is the signature to test for -- a column of identical values
    makes no ordering claim, so requiring a port to reproduce an ordering from
    it would be requiring it to agree about nothing.
    """
    return len({r["Ranking"] for r in rows.values()}) > 1


def _ranking_disagreements(golden, candidate, names, atol_es=ATOL_ES):
    """Count pairs the golden orders definitely but the candidate orders other.

    "Definitely" means the pair is distinguishable by the keys
    ``dpds_less_key`` sorts on: ``p_permutation`` first, then ``ES``
    descending. A pair that is not distinguishable has no defined order, so its
    ``Ranking`` cannot be required to match.

    The threshold for "not distinguishable" on ``ES`` is ``atol_es``, not exact
    equality, and that is the point. ``_numeric_deviation`` has already
    declared two scores within ``ATOL_ES`` to be the same number; demanding an
    ordering derived from their difference would contradict that in the same
    report. Found the hard way -- a faithful reimplementation agreed on every
    value and still failed two cases, because an ES difference of 2.776e-17,
    accepted as zero by the numeric check, swapped two rows and so counted as
    an ordering disagreement.

    ``p_permutation`` is compared exactly, because it is required to be exact.

    Returns ``(disagreements, n_tied_pairs)``.
    """
    disagreements = 0
    tied = 0
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            a, b = names[i], names[j]
            pa, pb = golden[a]["p_permutation"], golden[b]["p_permutation"]
            ea, eb = golden[a]["ES"], golden[b]["ES"]
            if pa == pb and abs(ea - eb) <= atol_es:
                tied += 1
                continue
            golden_first = (pa, -ea) < (pb, -eb)
            cand_first = candidate[a]["Ranking"] < candidate[b]["Ranking"]
            if golden_first != cand_first:
                disagreements += 1
    return disagreements, tied


def _hit_set(rows, threshold):
    return {n for n, r in rows.items() if r["FDR"] <= threshold}


def compare(case_name, golden_text, candidate_text, rtol=RTOL_PRINTED):
    """Compare two ``mageckGSEA`` outputs and return a :class:`Report`."""
    golden, golden_order = parse(golden_text)
    candidate, candidate_order = parse(candidate_text)

    report = Report(case=case_name)
    report.byte_identical = golden_text == candidate_text
    report.missing = sorted(set(golden) - set(candidate))
    report.extra = sorted(set(candidate) - set(golden))
    shared = sorted(set(golden) & set(candidate))

    for name in shared:
        for column in INTEGER_COLUMNS:
            g, c = golden[name][column], candidate[name][column]
            if g != c:
                report.structural.append(
                    "%s %s: golden=%r candidate=%r"
                    % (name or "<empty name>", column, g, c))
        # LFC is a hardcoded zero in both output branches. It is checked rather
        # than ignored so that a port inventing a value for it is caught.
        if candidate[name]["LFC"] != 0.0:
            report.structural.append(
                "%s LFC: candidate wrote %r, but the column is a literal 0"
                % (name or "<empty name>", candidate[name]["LFC"]))

    for column, printed, atol in NUMERIC_COLUMNS:
        report.deviations[column] = _numeric_deviation(
            column, golden, candidate, shared,
            rtol if printed else RTOL_PERMUTATION, atol)

    golden_ranked = is_ranked(golden)
    candidate_ranked = is_ranked(candidate)
    report.ranking_checked = golden_ranked and candidate_ranked
    if golden_ranked != candidate_ranked:
        report.structural.append(
            "Ranking: golden %s an ordering and candidate %s, so one of them "
            "took the sorted output branch and the other did not"
            % ("carries" if golden_ranked else "carries no",
               "does" if candidate_ranked else "does not"))

    if shared:
        gp = np.array([golden[n]["p_permutation"] for n in shared])
        cp = np.array([candidate[n]["p_permutation"] for n in shared])
        report.spearman = spearman(gp, cp)
        report.rank_inversions = inversions(gp, cp)
        if report.ranking_checked:
            report.ranking_disagreements, report.tied_pairs = (
                _ranking_disagreements(golden, candidate, shared))
    report.order_matches = golden_order == candidate_order

    for thr in HIT_THRESHOLDS:
        hg = _hit_set(golden, thr)
        hc = _hit_set(candidate, thr)
        report.hits[thr] = {
            "n_golden": len(hg),
            "n_candidate": len(hc),
            "only_golden": sorted(hg - hc),
            "only_candidate": sorted(hc - hg),
        }
    return report
