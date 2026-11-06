"""The fixed corpus of RRA inputs the port is validated against.

Every case here exists to cover a specific code path, named in its ``why``
field, and the corpus is chosen to cover the paths that *were not* covered by
the prototype benchmarking that preceded the port:

* The benchmark ran with ``--no-permutation-by-group``. That is not the
  default. ``heterogeneous`` and ``homogeneous_5`` exercise the by-group path,
  where a separate null is built for each distinct sgRNA-per-gene count.
* The benchmark used continuous values, so no two sgRNAs ever tied.
  ``bTreeSearchingF`` (``math_api.cpp``) resolves ties by recursion path rather
  than by any stated rule, which makes tie handling the single highest-risk
  part of the port. ``tie_heavy`` and ``tie_extreme`` target it directly.
* The benchmark never passed ``--control``. The control path replaces the
  uniform null with resampled control-sgRNA percentiles (``RRA.cpp:755-765``),
  which is a different null distribution, not a different way of computing the
  same one.

Inputs are generated from fixed seeds and then committed, rather than
regenerated at test time. NumPy's ``default_rng`` stream is stable by policy,
but committing the data removes the question entirely: a golden and the input
that produced it travel together.
"""

import dataclasses

import numpy as np

from .modes import Mode

#: Columns of the RRA input file, as ``fileio.cpp`` parses them:
#: ``<item id> <group id> <list id> <value>`` behind a header row.
INPUT_HEADER = "sgrna\tsymbol\tpool\tvalue"


@dataclasses.dataclass(frozen=True)
class Case:
    """One corpus entry: how to build it, how to run it, why it exists."""

    name: str
    why: str
    mode: Mode
    #: Extra RRA flags beyond ``-i``/``-o``/``-p``. ``--control`` is added by
    #: the runner when the case writes a control file, so it is not listed.
    flags: tuple = ()
    alpha: float = 0.25
    ngenes: int = 300
    #: An int for a uniform library, or a tuple cycled over genes to produce
    #: several permutation strata.
    sgrnas_per_gene: object = 5
    seed: int = 11
    hit_fraction: float = 0.02
    effect: float = -3.0
    #: Decimal places to round values to. Rounding is what manufactures exact
    #: ties; ``None`` leaves values continuous.
    round_to: object = None
    #: Fraction of genes whose sgRNAs are written to ``control.txt``.
    control_fraction: float = 0.0
    #: Write a control file whose names match nothing in the input. Used to
    #: reach the ``n_control == 0`` path.
    control_unmatched: bool = False
    #: Shape of the group column. ``"single"`` is one group name per row.
    #: ``"pathway"`` is the shape mageck2's pathway mode writes -- one row per
    #: gene, the group field a comma-separated list of every pathway the gene
    #: belongs to. ``"duplicate"`` repeats a name inside that list.
    group_style: str = "single"
    #: ``"pathway"``/``"duplicate"`` only: how many pathways to draw from, and
    #: the most any one gene joins.
    npathways: int = 40
    max_pathways_per_gene: int = 4
    #: ``"pathway"``/``"duplicate"`` only: how many pathways get a real effect.
    enriched_pathways: int = 3

    @property
    def uses_control(self):
        return self.control_fraction > 0.0 or self.control_unmatched


def _sgrna_counts(case):
    """Number of sgRNAs for each gene, as a list of length ``ngenes``."""
    spec = case.sgrnas_per_gene
    if isinstance(spec, int):
        return [spec] * case.ngenes
    return [spec[i % len(spec)] for i in range(case.ngenes)]


def _build_pathway(case, directory):
    """Write an input in the shape mageck2's pathway mode writes.

    ``pathwayFunc.py:146-162`` emits one row per gene, with the group field a
    comma-separated list of every pathway that contains the gene, and the
    literal ``NA`` for a gene in none. A gene in several pathways is therefore
    one row that belongs to several groups, which the C++ handles by splitting
    the field (``fileio.cpp:107`` and ``:266``) -- the group column is a list,
    not a name.

    This shape is the reason the case exists: it is what the shipped pathway
    mode actually produces, and nothing else in the corpus has a comma in the
    group column, so nothing else distinguishes "belongs to P1 and P2" from
    "belongs to a group called P1,P2".

    Returns the names of the pathways given a real effect.
    """
    rng = np.random.default_rng(case.seed)
    pathways = ["PATHWAY%03d" % i for i in range(case.npathways)]
    enriched = sorted(rng.choice(case.npathways, size=case.enriched_pathways,
                                 replace=False).tolist())
    enriched_names = {pathways[i] for i in enriched}

    rows = [INPUT_HEADER]
    for gene in range(case.ngenes):
        n_join = int(rng.integers(0, case.max_pathways_per_gene + 1))
        if n_join == 0:
            membership = []
        else:
            membership = sorted(rng.choice(
                case.npathways, size=n_join, replace=False).tolist())
        names = [pathways[i] for i in membership]
        if case.group_style == "duplicate" and len(names) >= 2:
            # A repeated sub-name is not collapsed: the C++ increments itemNum
            # once per occurrence, so the row lands in that group twice. Only
            # a corpus entry that repeats one proves the port agrees.
            names = [names[0]] + names
        # Genes in an enriched pathway get the effect; the rest are null.
        shift = case.effect if (enriched_names & set(names)) else 0.0
        value = rng.normal(shift, 1.0)
        if case.round_to is not None:
            value = round(value, case.round_to)
        field = ",".join(names) if names else "NA"
        rows.append("GENE%04d\t%s\tlist\t%.6f" % (gene, field, value))

    (directory / "input.txt").write_text("\n".join(rows) + "\n")
    return enriched_names


def build(case, directory):
    """Write ``input.txt`` (and ``control.txt`` if needed) for ``case``.

    Returns the set of group names given a real effect, so a test can check
    that an implementation recovers them without having to re-derive the
    planting. For the pathway-shaped cases those names are pathways, because
    pathways are what the output has one row per.
    """
    if case.group_style in ("pathway", "duplicate"):
        return _build_pathway(case, directory)
    rng = np.random.default_rng(case.seed)
    counts = _sgrna_counts(case)
    nhits = max(1, int(round(case.ngenes * case.hit_fraction)))
    hits = sorted(rng.choice(case.ngenes, size=nhits, replace=False).tolist())
    hitset = set(hits)

    names = ["GENE%04d" % g for g in range(case.ngenes)]
    rows = [INPUT_HEADER]
    control_names = []
    # Control genes are taken from the end of the list so they never collide
    # with the planted hits, which are drawn over the whole range.
    n_control_genes = int(round(case.ngenes * case.control_fraction))
    control_genes = set(range(case.ngenes - n_control_genes, case.ngenes)) - hitset

    for gene, k in enumerate(counts):
        shift = case.effect if gene in hitset else 0.0
        for i in range(k):
            value = rng.normal(shift, 1.0)
            if case.round_to is not None:
                value = round(value, case.round_to)
            sgrna = "sg%04d_%d" % (gene, i)
            rows.append("%s\t%s\tlist\t%.6f" % (sgrna, names[gene], value))
            if gene in control_genes:
                control_names.append(sgrna)

    (directory / "input.txt").write_text("\n".join(rows) + "\n")

    if case.control_unmatched:
        # Names that are syntactically plausible but absent from the input, so
        # every lookup in ProcessGroups misses and n_control ends at 0.
        control_names = ["absent_sgrna_%d" % i for i in range(20)]
    if case.uses_control:
        (directory / "control.txt").write_text("\n".join(control_names) + "\n")

    return {names[g] for g in hits}


#: Heterogeneous libraries in the wild are not evenly spread over sgRNA counts;
#: most genes carry the library's nominal number and a tail carries fewer
#: because of filtering. This mix produces eight permutation strata with
#: uneven occupancy, which is the shape that actually occurs.
HETEROGENEOUS_MIX = (4, 4, 4, 4, 5, 5, 5, 6, 3, 2, 1, 8, 12, 4, 5, 4)

CASES = (
    Case(
        name="homogeneous_5",
        why="Baseline. One permutation stratum, the default by-group path, "
            "continuous values so no ties. If this case fails nothing else "
            "is worth reading.",
        mode=Mode.PARITY,
        sgrnas_per_gene=5,
    ),
    Case(
        name="heterogeneous",
        why="The default path on a realistic library: eight distinct "
            "sgRNA-per-gene counts means eight separate nulls, each sized by "
            "its own stratum. The prototype benchmark never exercised this.",
        mode=Mode.PARITY,
        sgrnas_per_gene=HETEROGENEOUS_MIX,
        seed=23,
    ),
    Case(
        name="tie_heavy",
        why="Values rounded to two decimals, so sgRNAs tie by the dozen. "
            "Percentiles come from bTreeSearchingF, whose tie resolution "
            "depends on the recursion path rather than on a stated rule, so "
            "this is where a port most plausibly diverges.",
        mode=Mode.PARITY,
        sgrnas_per_gene=HETEROGENEOUS_MIX,
        round_to=2,
        seed=7,
    ),
    Case(
        name="tie_extreme",
        why="Values rounded to whole numbers, which collapses a 300-gene "
            "screen onto about a dozen distinct values. Degenerate on "
            "purpose: it turns a rare tie-breaking disagreement into one that "
            "shows up in almost every gene.",
        mode=Mode.PARITY,
        sgrnas_per_gene=5,
        round_to=0,
        seed=31,
    ),
    Case(
        name="pooled_permutation",
        why="--no-permutation-by-group: one null for all genes regardless of "
            "sgRNA count. This is the path the prototype benchmark validated, "
            "kept so a regression in it is still caught.",
        mode=Mode.PARITY,
        flags=("--no-permutation-by-group",),
        sgrnas_per_gene=HETEROGENEOUS_MIX,
        seed=23,
    ),
    Case(
        name="control_matched",
        why="--control with a list that matches real sgRNAs. The null is "
            "resampled from control-sgRNA percentiles instead of drawn "
            "uniformly, which is a different distribution rather than a "
            "different route to the same one. Never previously validated.",
        mode=Mode.PARITY,
        sgrnas_per_gene=5,
        control_fraction=0.15,
        seed=41,
    ),
    Case(
        name="control_heterogeneous",
        why="The control null combined with by-group permutation, so the "
            "resampled percentiles are drawn once per sgRNA count. The two "
            "features interact and are worth testing together.",
        mode=Mode.PARITY,
        sgrnas_per_gene=HETEROGENEOUS_MIX,
        control_fraction=0.15,
        seed=43,
    ),
    Case(
        name="single_sgrna_genes",
        why="Many genes with exactly one sgRNA. ComputeLoValue always "
            "evaluates index 0 before applying the alpha cutoff, so a "
            "one-sgRNA gene takes a special path through the kept-set logic.",
        mode=Mode.PARITY,
        sgrnas_per_gene=(1, 1, 1, 2, 5),
        seed=53,
    ),
    Case(
        name="alpha_low",
        why="-p 0.05, the bottom of the range crisprFunction.py clamps to. A "
            "small alpha keeps few sgRNAs per gene, so goodsgrna counts and "
            "the kept-set prefix logic are stressed differently.",
        mode=Mode.PARITY,
        alpha=0.05,
        sgrnas_per_gene=HETEROGENEOUS_MIX,
        seed=23,
    ),
    Case(
        name="alpha_high",
        why="-p 0.5, the top of that range, where nearly every sgRNA is kept.",
        mode=Mode.PARITY,
        alpha=0.5,
        sgrnas_per_gene=HETEROGENEOUS_MIX,
        seed=23,
    ),
    Case(
        name="min_goodsgrna",
        why="--min-number-goodsgrna excludes genes from FDR adjustment "
            "(adjustFDR, RRA.cpp:646), which makes them print FDR as the "
            "string NA. The port has to reproduce which rows are NA, not just "
            "the numbers in the rows that are not.",
        mode=Mode.PARITY,
        flags=("--min-number-goodsgrna", "3"),
        sgrnas_per_gene=HETEROGENEOUS_MIX,
        seed=23,
    ),
    Case(
        name="skip_gene",
        why="--skip-gene, which crisprFunction.py always passes (as NA/na). "
            "Skipped genes bypass permutation, so they exercise a separate "
            "branch of the group loop.",
        mode=Mode.PARITY,
        flags=("--skip-gene", "GENE0000", "--skip-gene", "GENE0001"),
        sgrnas_per_gene=5,
    ),
    Case(
        name="permutation_ceiling",
        why="--max-permutation-per-group bounds the null size and so raises "
            "the smallest reachable p value. The relationship between that "
            "bound and p resolution is asserted in test_rra.py; this case "
            "checks the port reproduces the bounded null itself.",
        mode=Mode.PARITY,
        flags=("--max-permutation-per-group", "20000"),
        sgrnas_per_gene=HETEROGENEOUS_MIX,
        seed=23,
    ),
    Case(
        name="min_percentile",
        why="-P, the lower percentile bound. Nothing in the corpus passed it "
            "before, because its default of -1.0 cannot exclude anything: the "
            "port's lower bound was therefore unreachable, and every golden "
            "recorded so far is consistent with it being wrong. -P 0.01 "
            "excludes the strongest sgRNAs -- except at index 0, which "
            "ComputeLoValue scores whether or not it is in range "
            "(RRA.cpp:500-502), so a gene can end up with a lo-value below 1 "
            "and no good sgRNAs at all. This is the case that fails if the "
            "bound is applied to the observed lo-values but not to the "
            "permutation null: 25 genes move outside the FDR tolerance.",
        mode=Mode.PARITY,
        flags=("-P", "0.01"),
        sgrnas_per_gene=HETEROGENEOUS_MIX,
        seed=23,
    ),
    Case(
        name="min_percentile_control",
        why="-P with --control, so the bound applies to a null resampled from "
            "control-sgRNA percentiles rather than a uniform one. The C++ "
            "passes minPercentile into ComputePermutationbyNumSG as well as "
            "into ProcessGroups (RRA.cpp:306-308), so the bound has to reach "
            "the null too.\n\n"
            "The bound is 0.08 rather than 0.01 because it was chosen from "
            "the data instead of by eye. This case's smallest control "
            "percentile is 0.0215, so -P 0.01 excludes nothing a resampled "
            "null can draw: the case passed whether or not the null was "
            "bounded, which makes it no test of the thing it is named for. "
            "0.08 excludes 14 of the 206 control percentiles, so an unbounded "
            "null now shows up as wrong p values.",
        mode=Mode.PARITY,
        flags=("-P", "0.08"),
        sgrnas_per_gene=HETEROGENEOUS_MIX,
        control_fraction=0.15,
        seed=41,
    ),
    Case(
        name="pathway_memberships",
        why="The group column as the shipped pathway mode writes it: one row "
            "per gene, the field a comma-separated list of every pathway "
            "containing that gene, NA for none (pathwayFunc.py:155-162). The "
            "C++ splits the field and counts the row in each pathway. Nothing "
            "else in the corpus has a comma in that column, which is how the "
            "first port of this module came to treat P1,P2 as a group named "
            "\"P1,P2\" -- inventing a pathway and leaving the two real ones "
            "a member short each -- and why that went unnoticed.",
        mode=Mode.PARITY,
        group_style="pathway",
        ngenes=600,
        npathways=40,
        seed=101,
    ),
    Case(
        name="pathway_duplicate_membership",
        why="Same shape, but with a sub-name repeated inside the field. The "
            "C++ does not deduplicate: it appends the row once per "
            "occurrence, so P1,P1 gives P1 the row twice and its sgRNA count "
            "is larger than the number of distinct rows in it. Measured "
            "against the helper, not assumed.",
        mode=Mode.PARITY,
        group_style="duplicate",
        ngenes=400,
        npathways=25,
        seed=103,
    ),
    Case(
        name="control_empty",
        why="--control with a list matching nothing. n_control ends at 0, so "
            "RRA.cpp:761-762 computes rand_ctl_index = -1 and reads "
            "control_prob_array[-1] -- the out-of-bounds read in issue #51. "
            "The C++ result here is undefined, so there is nothing to copy: "
            "the port must choose a behaviour and this case pins the choice.",
        mode=Mode.CONTRACT,
        sgrnas_per_gene=5,
        control_unmatched=True,
    ),
)

BY_NAME = {c.name: c for c in CASES}

PARITY_CASES = tuple(c for c in CASES if c.mode is Mode.PARITY)
CONTRACT_CASES = tuple(c for c in CASES if c.mode is Mode.CONTRACT)
