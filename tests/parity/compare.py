"""Column-aware comparison of two RRA outputs.

Comparing these files byte for byte would be the wrong test in both
directions. It is too strict, because ``fileio.cpp:337`` prints ``lo_value``
and ``p`` through ``%10.4e`` -- five significant figures -- so the file already
discards precision the implementations may legitimately differ in below. And it
is too weak, because a byte-identical diff says nothing about *which* genes
moved or whether anyone would notice: rows are written in ``lo_value`` order,
so a single swapped pair of tied genes shifts every subsequent line and the
diff looks catastrophic when nothing scientifically changed.

So this module compares by gene, applies a rule per column, and reports three
levels of agreement:

1. **Structural** -- the same genes are present, with the same sgRNA counts and
   the same set of rows excluded from FDR adjustment.
2. **Numeric** -- each column agrees within a tolerance matched to the
   precision the file actually carries, with the worst observed deviation
   reported even when it passes, so a drift that stays inside tolerance is
   still visible.
3. **Scientific** -- the p value ranking and the called hit sets agree. This is
   the level a user would notice, and it is reported separately because a port
   can pass (2) and still fail (3) if small deviations straddle an FDR
   threshold where many genes are packed together.
"""

import dataclasses

import numpy as np

from .stats import inversions as _inversions
from .stats import normalized_ranks  # noqa: F401  (re-exported for tests)
from .stats import spearman as _spearman

#: Columns of the RRA output, in the order ``fileio.cpp:329`` writes them.
COLUMNS = ("group_id", "items_in_group", "lo_value", "p", "FDR", "goodsgrna")

#: ``%10.4e`` carries five significant figures, so two implementations that
#: agree to the last printed digit differ by at most ~5e-5 relatively. 1e-4
#: accepts that and nothing more.
RTOL_PRINTED = 1e-4

#: ``FDR`` is printed through ``%f``: six decimals, fixed point. The tolerance
#: is a little looser than half an ulp because FDR comes out of a step-down
#: pass over ranked p values, so a p difference far below tolerance can still
#: move an FDR value by a whole step when p values are tied or near-tied.
ATOL_FDR = 5e-6

#: FDR cutoffs at which the called hit sets are compared. 0.05 and 0.1 are the
#: conventional reporting thresholds; 0.25 is the MAGeCK default alpha and the
#: value `--gene-test-fdr-threshold` defaults to.
HIT_THRESHOLDS = (0.05, 0.1, 0.25)


def parse(text):
    """Parse an RRA output file into ``{gene: {column: value}}``.

    ``lo_value``/``p`` become floats, ``items_in_group``/``goodsgrna`` ints,
    and ``FDR`` either a float or ``None`` for the literal ``NA`` that
    ``fileio.cpp:339`` writes for genes held out of FDR adjustment.
    """
    lines = [ln for ln in text.splitlines() if ln.strip()]
    if not lines:
        raise ValueError("empty RRA output")
    header = [f.strip() for f in lines[0].split("\t")]
    if tuple(header) != COLUMNS:
        raise ValueError("unexpected RRA output header: %r" % (header,))

    rows = {}
    order = []
    for ln in lines[1:]:
        f = [x.strip() for x in ln.split("\t")]
        if len(f) != len(COLUMNS):
            raise ValueError("malformed RRA output row: %r" % (ln,))
        gene = f[0]
        rows[gene] = {
            "items_in_group": int(f[1]),
            "lo_value": float(f[2]),
            "p": float(f[3]),
            "FDR": None if f[4] == "NA" else float(f[4]),
            "goodsgrna": int(f[5]),
        }
        order.append(gene)
    return rows, order


@dataclasses.dataclass
class Deviation:
    """Worst observed disagreement in one numeric column."""

    column: str
    max_abs: float = 0.0
    max_rel: float = 0.0
    worst_gene: str = ""
    n_outside: int = 0

    def __str__(self):
        s = "%-14s max_abs=%.3e max_rel=%.3e" % (
            self.column, self.max_abs, self.max_rel)
        if self.worst_gene:
            s += "  worst=%s" % self.worst_gene
        if self.n_outside:
            s += "  OUTSIDE TOLERANCE: %d gene(s)" % self.n_outside
        return s


@dataclasses.dataclass
class Report:
    """The outcome of comparing a candidate against a golden."""

    case: str
    missing: list = dataclasses.field(default_factory=list)
    extra: list = dataclasses.field(default_factory=list)
    structural: list = dataclasses.field(default_factory=list)
    na_only_golden: list = dataclasses.field(default_factory=list)
    na_only_candidate: list = dataclasses.field(default_factory=list)
    deviations: dict = dataclasses.field(default_factory=dict)
    spearman: float = float("nan")
    rank_inversions: int = 0
    order_matches: bool = True
    hits: dict = dataclasses.field(default_factory=dict)

    @property
    def ok(self):
        """True when every asserted level agrees.

        Row order is deliberately excluded: it is reported but not required,
        because rows are sorted on a value that legitimately ties.
        """
        return (not self.missing and not self.extra and not self.structural
                and not self.na_only_golden and not self.na_only_candidate
                and all(d.n_outside == 0 for d in self.deviations.values())
                and all(h["only_golden"] == [] and h["only_candidate"] == []
                        for h in self.hits.values()))

    def render(self):
        out = ["case %s: %s" % (self.case, "AGREE" if self.ok else "DISAGREE")]
        if self.missing:
            out.append("  genes missing from candidate: %d %s"
                       % (len(self.missing), self.missing[:8]))
        if self.extra:
            out.append("  genes only in candidate: %d %s"
                       % (len(self.extra), self.extra[:8]))
        for line in self.structural[:12]:
            out.append("  " + line)
        if len(self.structural) > 12:
            out.append("  ... %d more structural mismatches"
                       % (len(self.structural) - 12))
        if self.na_only_golden:
            out.append("  FDR is NA only in golden: %d %s"
                       % (len(self.na_only_golden), self.na_only_golden[:8]))
        if self.na_only_candidate:
            out.append("  FDR is NA only in candidate: %d %s"
                       % (len(self.na_only_candidate),
                          self.na_only_candidate[:8]))
        for dev in self.deviations.values():
            out.append("  " + str(dev))
        out.append("  p ranking: spearman=%.6f inversions=%d row_order=%s"
                   % (self.spearman, self.rank_inversions,
                      "same" if self.order_matches else "differs"))
        for thr, h in sorted(self.hits.items()):
            out.append("  hits at FDR<=%.2f: golden=%d candidate=%d "
                       "only_golden=%s only_candidate=%s"
                       % (thr, h["n_golden"], h["n_candidate"],
                          h["only_golden"][:6], h["only_candidate"][:6]))
        return "\n".join(out)


def _numeric_deviation(column, golden, candidate, genes, rtol, atol):
    dev = Deviation(column=column)
    for gene in genes:
        g = golden[gene][column]
        c = candidate[gene][column]
        if g is None or c is None:
            # NA handling is structural, reported separately.
            continue
        abs_err = abs(g - c)
        if g != 0.0:
            rel_err = abs_err / abs(g)
        else:
            rel_err = 0.0 if c == 0.0 else float("inf")
        if abs_err > dev.max_abs:
            dev.max_abs = abs_err
            dev.worst_gene = gene
        dev.max_rel = max(dev.max_rel, rel_err)
        if abs_err > atol + rtol * abs(g):
            dev.n_outside += 1
    return dev


def _hit_set(rows, threshold):
    return {g for g, r in rows.items()
            if r["FDR"] is not None and r["FDR"] <= threshold}


def compare(case_name, golden_text, candidate_text,
            rtol=RTOL_PRINTED, atol_fdr=ATOL_FDR):
    """Compare two RRA outputs and return a :class:`Report`."""
    golden, golden_order = parse(golden_text)
    candidate, candidate_order = parse(candidate_text)

    report = Report(case=case_name)
    report.missing = sorted(set(golden) - set(candidate))
    report.extra = sorted(set(candidate) - set(golden))
    shared = sorted(set(golden) & set(candidate))

    for gene in shared:
        for column in ("items_in_group", "goodsgrna"):
            g, c = golden[gene][column], candidate[gene][column]
            if g != c:
                report.structural.append(
                    "%s %s: golden=%r candidate=%r" % (gene, column, g, c))
        g_na = golden[gene]["FDR"] is None
        c_na = candidate[gene]["FDR"] is None
        if g_na and not c_na:
            report.na_only_golden.append(gene)
        elif c_na and not g_na:
            report.na_only_candidate.append(gene)

    for column, atol in (("lo_value", 0.0), ("p", 0.0), ("FDR", atol_fdr)):
        report.deviations[column] = _numeric_deviation(
            column, golden, candidate, shared, rtol, atol)

    if shared:
        gp = np.array([golden[g]["p"] for g in shared])
        cp = np.array([candidate[g]["p"] for g in shared])
        report.spearman = _spearman(gp, cp)
        report.rank_inversions = _inversions(gp, cp)
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
