"""The fixed corpus of ``mageckGSEA`` inputs the port is validated against.

Each case names the code path it exists to cover in its ``why`` field. The
corpus is built around three things the existing ``tests/test_gsea.py`` does
not reach:

* **Both rank-file branches.** ``mageck pathway`` always passes ``-c 2`` or
  ``-c 8``, so the ``-c 0`` branch -- where a gene's score comes from its line
  number rather than from a column -- has no coverage at all even though it is
  a documented option and a completely separate scoring rule.

* **The tie regime.** With ``-p 1000`` and a few dozen pathways, permutation p
  values land on a 1/1000 grid and pathways tie by the handful. Ties decide
  both the output row order and, through the step-down pass, the FDR. The
  ``tie_*`` and ``coarse_permutation`` cases make that the common case rather
  than an accident.

* **Degenerate gene sets.** A set covering the whole ranked list, a set
  overlapping it not at all, a set of duplicates, and a set of one each take a
  separate early return in ``getgseascore_fast`` or ``getscoreandp``. Real GMT
  files contain all four.

Inputs are generated from fixed seeds and then committed, rather than
regenerated at test time, so that a golden and the input that produced it
travel together.
"""

import dataclasses

import numpy as np

from ..modes import Mode

#: Header of the generated rank file. Shaped like the ``gene_summary.txt`` that
#: ``mageck pathway`` is actually pointed at, because the column index passed
#: with ``-c`` only means anything relative to a real layout: 2 is the negative
#: selection score and 8 the positive one (``argsParser.py``, ``--ranking-column``
#: and ``--ranking-column-2``).
RANK_HEADER = ("id\tnum\tneg|score\tneg|p-value\tneg|fdr\tneg|rank\t"
               "neg|goodsgrna\tneg|lfc\tpos|score\tpos|p-value")

#: Number of columns ``RANK_HEADER`` declares. Asserted against the rows
#: :func:`build` writes, because those are assembled from a hand-written list
#: of ten fields and a header that drifted out of step with them would move
#: every column index -- including the ``-c 2`` and ``-c 8`` the production
#: call uses -- while still producing a file that parses.
RANK_COLUMNS = len(RANK_HEADER.split("\t"))

#: The GMT description column written for every generated gene set. One token,
#: for the reason given in :func:`build`.
DESCRIPTION = "generated_gene_set"


@dataclasses.dataclass(frozen=True)
class Case:
    """One corpus entry: how to build it, how to run it, why it exists."""

    name: str
    why: str
    mode: Mode

    # --- how the helper is invoked ---------------------------------------
    #: ``-c``: which rank-file column scores the genes. 0 selects the
    #: line-number branch of ``parserankfile``.
    score_column: int = 2
    #: ``-p``: permutations. 1000 is what ``mageck pathway`` passes.
    perms: int = 200
    #: ``--seed``. ``None`` leaves it unset, which is the production call and
    #: exercises the default.
    seed: object = None
    #: ``-s``: sort by p value. ``mageck pathway`` always passes it; without it
    #: the output takes a different branch that writes p into the FDR column.
    sort_by_p: bool = True
    #: ``-H``: skip the rank file's header row. Always passed in production.
    skip_header: bool = True
    #: ``-e``: reverse the gene order after scoring -- and only the names, not
    #: the scores, which is why it is not a symmetry.
    reverse: bool = False
    #: ``-n``: test one named pathway instead of all of them. Takes a separate
    #: branch in ``main``.
    pathway_name: object = None

    # --- how the input is generated --------------------------------------
    ngenes: int = 400
    nsets: int = 30
    #: Gene-set sizes are drawn uniformly from this inclusive range.
    set_size: tuple = (5, 60)
    #: Every ``enrich_every``-th set is drawn from the top of the ranking
    #: instead of uniformly, so some sets are genuinely enriched and the
    #: corpus is not all null.
    enrich_every: int = 5
    #: Decimal places to round scores to. Rounding is what manufactures the
    #: exact ties that decide sort order; ``None`` leaves them continuous.
    round_to: object = 6
    gen_seed: int = 11

    # --- deliberate input malformations ----------------------------------
    #: Omit the final newline of the GMT file.
    gmt_no_trailing_newline: bool = False
    #: Omit the final newline of the rank file.
    rank_no_trailing_newline: bool = False
    #: Insert a blank line into the middle of the GMT file.
    gmt_blank_line: bool = False
    #: Leave trailing spaces on the first GMT line.
    gmt_trailing_whitespace: bool = False
    #: Give the first gene set a description containing spaces, one of whose
    #: words is a real gene name from the top of the ranking that is *not* a
    #: member of that set.
    gmt_description_with_spaces: bool = False
    #: Add a set containing every gene in the ranked list.
    whole_list_set: bool = False
    #: Add a set whose genes are all absent from the ranked list.
    absent_set: bool = False
    #: Add a set that names the same gene several times.
    duplicate_set: bool = False
    #: Add several sets of exactly one gene.
    singleton_sets: bool = False

    #: Whether the enriched sets this case plants are expected to rank above
    #: the null sets in *this* case's output. They are not always: the planting
    #: is relative to the ``neg|score`` ordering, so scoring on ``pos|score``
    #: (``-c 8``) puts them at the other end, and the ``-c 0`` branch ignores
    #: the score columns altogether. Stated per case rather than inferred,
    #: because a test that guesses this would be asserting the generator's
    #: behaviour rather than the port's.
    signal_expected: bool = True

    def gene_names(self):
        return ["GENE%05d" % i for i in range(self.ngenes)]


def _set_members(case, rng, genes):
    """The gene sets, as an ordered list of ``(name, members)``.

    Built in a fixed order; :func:`build` then shuffles that order before
    writing, so that file order and the lexicographic order the helper
    actually iterates in are different things. See the comment there.
    """
    sets = []
    lo, hi = case.set_size
    for i in range(case.nsets):
        size = int(rng.integers(lo, hi + 1))
        if case.enrich_every and i % case.enrich_every == 0:
            pool = genes[:max(size, len(genes) // 8)]
        else:
            pool = genes
        size = min(size, len(pool))
        members = rng.choice(len(pool), size=size, replace=False)
        sets.append(("SET_%03d" % i, [pool[j] for j in sorted(members)]))

    if case.whole_list_set:
        # nf == n in getgseascore_fast, so decvalue would divide by zero.
        sets.append(("SET_WHOLE_LIST", list(genes)))
    if case.absent_set:
        # Every GNAMEINDEX lookup misses, so index stays empty.
        sets.append(("SET_ABSENT", ["NOT_IN_RANKING_%d" % i for i in range(6)]))
    if case.duplicate_set:
        # checkred dedups, so Size (the raw GMT count) and the number of genes
        # actually scored come apart.
        sets.append(("SET_DUPLICATED",
                     [genes[3], genes[3], genes[3], genes[7], genes[7]]))
    if case.singleton_sets:
        # num == 1, so gindex is always 0 and p_gsea is betainc(1, 1, x) == x.
        for i, g in enumerate((genes[0], genes[len(genes) // 2], genes[-1])):
            sets.append(("SET_SINGLE_%d" % i, [g]))
    return sets


def build(case, directory):
    """Write ``sets.gmt`` and ``rank.txt`` for ``case``.

    Returns the names of the gene sets that were drawn from the top of the
    ranking, so a test can check that an implementation ranks them above the
    null sets without having to re-derive the planting.
    """
    rng = np.random.default_rng(case.gen_seed)
    genes = case.gene_names()

    # Scores are drawn uniformly, so the ranking is non-trivial. The file is
    # not written in score order either -- the helper sorts it itself, and a
    # port that relied on the input being sorted has to fail.
    scores = rng.uniform(0.0, 1.0, size=case.ngenes)
    if case.round_to is not None:
        scores = np.round(scores, case.round_to)
    order = np.argsort(scores, kind="stable")
    ranked = [genes[i] for i in order]

    # Rows are written in a shuffled order, which matters more than it looks.
    # Written in gene-name order, the file would have file order == name order,
    # and parserankfile's tie-break -- it sorts (score, name) pairs -- would be
    # indistinguishable from a stable sort by score alone. A port that dropped
    # the name tie-break entirely would then agree with every golden. Found by
    # injecting exactly that fault into a faithful reimplementation and
    # watching the harness accept it.
    layout = list(rng.permutation(case.ngenes))

    rows = [RANK_HEADER]
    for line_number, i in enumerate(layout, start=1):
        # Columns beyond the two score columns are filled with plausible but
        # unused values, so that a port reading the wrong column produces
        # something visibly wrong rather than something subtly wrong.
        neg = scores[i]
        pos = 1.0 - neg
        rows.append("\t".join([
            genes[i], "5", "%.6f" % neg, "%.6f" % neg,
            "%.6f" % min(1.0, neg * 2), str(line_number), "4",
            "%.4f" % (neg - 0.5), "%.6f" % pos, "%.6f" % pos,
        ]))
    for row in rows:
        assert len(row.split("\t")) == RANK_COLUMNS, (
            "rank row has %d fields but RANK_HEADER declares %d: %r"
            % (len(row.split("\t")), RANK_COLUMNS, row[:60]))

    rank_text = "\n".join(rows)
    rank_text += "" if case.rank_no_trailing_newline else "\n"
    (directory / "rank.txt").write_text(rank_text)

    sets = _set_members(case, rng, ranked)

    # And the GMT lines are written in a shuffled order, for the same kind of
    # reason. The helper iterates a std::map, so it reads gene sets back in
    # lexicographic name order whatever order the file is in -- but with the
    # file already in name order, "lexicographic" and "file order" are the
    # same thing, and a port iterating the file would consume the permutation
    # stream identically and agree with every golden. Also found by injecting
    # the fault and watching the harness accept it.
    written = [sets[i] for i in rng.permutation(len(sets))]

    # A top-ranked gene that is not a member of the first *written* set, so
    # that reading it out of the description column visibly changes that set's
    # score rather than being a no-op.
    first_members = set(written[0][1]) if written else set()
    injected_gene = next((g for g in ranked if g not in first_members),
                         ranked[0])

    lines = []
    for i, (name, members) in enumerate(written):
        # A single-token description on purpose. GMT is a tab-delimited format
        # and its second column is free text, but parsegmt tokenises on
        # whitespace, so a description containing spaces has its second and
        # later words read as gene names. That is what the
        # gmt_description_with_spaces case is for; letting it happen in every
        # case would pollute every Size column and hide the finding.
        description = DESCRIPTION
        if case.gmt_description_with_spaces and i == 0:
            description = "response to %s" % injected_gene
        line = "\t".join([name, description] + list(members))
        if case.gmt_trailing_whitespace and i == 0:
            line += "   "
        lines.append(line)
    if case.gmt_blank_line and len(lines) > 1:
        lines.insert(1, "")
    gmt_text = "\n".join(lines)
    gmt_text += "" if case.gmt_no_trailing_newline else "\n"
    (directory / "sets.gmt").write_text(gmt_text)

    enriched = {name for i, (name, _) in enumerate(sets)
                if name.startswith("SET_") and name[4:].isdigit()
                and case.enrich_every and i % case.enrich_every == 0}
    return enriched


CASES = (
    Case(
        name="production_default",
        why="Exactly what `mageck pathway` runs for negative selection: "
            "-s -H -c 2 -p 1000. Every other case is a variation on this one, "
            "so if this fails nothing else is worth reading.",
        mode=Mode.PARITY,
        score_column=2,
        perms=1000,
    ),
    Case(
        name="production_positive",
        why="The second call `mageck pathway` makes, on the positive-selection "
            "column (-c 8). Same code path as the first, but it proves the "
            "column index is honoured rather than hardcoded -- a port reading "
            "a fixed column would pass production_default and fail here.",
        mode=Mode.PARITY,
        score_column=8,
        perms=1000,
        signal_expected=False,
    ),
    Case(
        name="rank_only",
        why="-c 0: the other branch of parserankfile entirely. The score comes "
            "from a gene's line number, the genes are not re-sorted, and the "
            "normalisation counts from 1 instead of 0. No existing test "
            "reaches this branch.",
        mode=Mode.PARITY,
        score_column=0,
        signal_expected=False,
    ),
    Case(
        name="unsorted",
        why="Without -s, testallgseas takes a different output branch that "
            "writes the permutation p into the FDR column and leaves Ranking "
            "at 0. Not what `mageck pathway` passes, but it is the documented "
            "default of the helper's own CLI.",
        mode=Mode.PARITY,
        sort_by_p=False,
    ),
    Case(
        name="no_header",
        why="Without -H the header row is parsed as a gene. atof() of a column "
            "label is 0.0, so a phantom gene joins the ranking and shifts the "
            "line count every other score is divided by. Reproducing this "
            "matters because -H is a flag a user can forget.",
        mode=Mode.PARITY,
        skip_header=False,
        signal_expected=False,
    ),
    Case(
        name="reverse",
        why="-e reverses the gene *names* after scoring but leaves the score "
            "vector alone (calculategsea.cpp:237-239), so it is not the "
            "symmetry its name suggests: gene i keeps position i's score. A "
            "port that reverses both agrees on nothing here.",
        mode=Mode.PARITY,
        reverse=True,
        signal_expected=False,
    ),
    Case(
        name="single_pathway",
        why="-n names one pathway, which takes the branch in main() rather "
            "than testallgseas, prints Ranking as 1 and runs no FDR pass at "
            "all. Also the only path where the permutation stream is consumed "
            "by one pathway instead of all of them in name order.",
        mode=Mode.PARITY,
        # A null set, not one of the enriched ones. An enriched set's
        # permutation p saturates at 0 -- no permutation ever beats it -- and a
        # saturated p cannot detect a change to the generator, so the case
        # would have covered the branch without covering the stream. Measured:
        # with SET_005 here, this was the one case in the corpus that a wrong
        # RNG multiplier did not move.
        pathway_name="SET_006",
        signal_expected=False,
    ),
    Case(
        name="tie_heavy",
        why="Scores rounded to two decimals, so genes tie by the dozen. Ties "
            "are resolved inside parserankfile by sorting (score, name) "
            "pairs, which means the gene order -- and therefore every "
            "score -- depends on a name comparison a port could easily omit.",
        mode=Mode.PARITY,
        round_to=2,
        gen_seed=23,
    ),
    Case(
        name="tie_extreme",
        why="Scores rounded to one decimal, which collapses 400 genes onto "
            "about ten distinct values. Degenerate on purpose: it turns a rare "
            "tie-breaking disagreement into one that shows up everywhere.",
        mode=Mode.PARITY,
        round_to=1,
        gen_seed=29,
    ),
    Case(
        name="coarse_permutation",
        why="-p 20, so permutation p values land on a 1/20 grid: 30 pathways "
            "share 14 distinct p values across 37 tied pairs, which the FDR "
            "step-down pass then collapses to 7 distinct FDR values. That is "
            "where a port's FDR can differ while every p agrees, because the "
            "step-down takes a running minimum over tied runs rather than "
            "adjusting each value independently. ES still breaks almost every "
            "tie, so the row ordering here stays defined -- the cases with "
            "genuinely unordered pairs are tie_extreme and big_sets.",
        mode=Mode.PARITY,
        perms=20,
    ),
    Case(
        name="degenerate_sets",
        why="A set covering the whole ranked list (n <= nf, the division "
            "getgseascore_fast guards), a set overlapping it not at all "
            "(empty index), a set of duplicates (Size and the scored count "
            "come apart), and three singletons (num == 1, so p_gsea reduces "
            "to the percentile itself). Real GMT files contain all four.",
        mode=Mode.PARITY,
        whole_list_set=True,
        absent_set=True,
        duplicate_set=True,
        singleton_sets=True,
    ),
    Case(
        name="tiny_sets",
        why="Every set has between one and three genes, so num - gindex is "
            "always small and BetaNoncentralCdf is evaluated at the edge of "
            "its parameter range rather than in the comfortable middle.",
        mode=Mode.PARITY,
        set_size=(1, 3),
        nsets=40,
        gen_seed=31,
    ),
    Case(
        name="big_sets",
        why="Sets covering a quarter to three quarters of the ranked list. The "
            "permutation draws then reject often -- sampling 300 distinct "
            "indices out of 400 by rejection takes far more than 300 draws -- "
            "so the RNG consumption pattern is completely different from the "
            "small-set case, and a port that samples without reproducing the "
            "rejection loop diverges only here.",
        mode=Mode.PARITY,
        set_size=(100, 300),
        nsets=12,
        gen_seed=37,
    ),
    Case(
        name="custom_seed",
        why="--seed 20260311. The seed reached the helper only in PR #52; "
            "before that it reseeded from the clock. This case fails if a port "
            "accepts the argument and ignores it, which is the quiet way to "
            "get that wrong.",
        mode=Mode.PARITY,
        seed=20260311,
    ),
    Case(
        name="many_pathways",
        why="150 sets over 2000 genes at -p 1000: the scale at which the "
            "permutation stream is consumed millions of times and in strict "
            "pathway-name order. A port that parallelises over pathways, or "
            "iterates them in file order rather than lexicographic order, "
            "produces plausible numbers that agree with nothing.",
        mode=Mode.PARITY,
        ngenes=2000,
        nsets=150,
        perms=1000,
        set_size=(5, 120),
        gen_seed=43,
    ),

    # --- malformed inputs: both implementations must be unmoved ----------
    Case(
        name="gmt_missing_final_newline",
        why="parsegmt loops `while(true){getline; if(eof) break;}`. On a file "
            "whose last line has no newline, that getline sets eofbit while "
            "returning the line, so the loop breaks before storing it: the "
            "last gene set is silently dropped. Confirmed against the binary "
            "-- it reports one pathway fewer and says nothing. The port must "
            "keep the set.",
        mode=Mode.PARITY,
        gmt_no_trailing_newline=True,
        nsets=8,
    ),
    Case(
        name="rank_missing_final_newline",
        why="The same loop in parserankfile, and worse: dropping the last gene "
            "also lowers the line count that every gene's score is divided by, "
            "so one missing newline perturbs the whole ranking rather than "
            "just losing a row. The port must keep the gene.",
        mode=Mode.PARITY,
        rank_no_trailing_newline=True,
        nsets=8,
    ),
    Case(
        name="gmt_blank_line",
        why="A blank line in the GMT makes parsegmt store a pathway whose name "
            "is the empty string. main() then finds PM.count(args.pname) > 0 "
            "-- because args.pname defaults to \"\" -- and takes the "
            "single-pathway branch, so it reports that one empty pathway and "
            "silently skips every real one, exiting 0. Confirmed against the "
            "binary: a 2-set file plus one blank line produces a single junk "
            "row. The port must ignore blank lines.",
        mode=Mode.PARITY,
        gmt_blank_line=True,
        nsets=8,
    ),
    Case(
        name="gmt_description_with_spaces",
        why="GMT is a tab-delimited format whose second column is free-text "
            "description, but parsegmt reads the line with `ss >> ssf` and "
            "treats every token after the second as a gene name. So a "
            "description containing spaces has its second and later words "
            "silently added to the gene set. Mostly they match nothing and "
            "only inflate Size -- but gene symbols occur in prose "
            "descriptions constantly (\"regulation of TP53 activity\"), and "
            "a word that does match is added to the set and changes its "
            "result. Confirmed against the binary: giving a two-gene set the "
            "description \"response to G01\" moved ES from 5.55e-17 to "
            "0.396226 and p from 0.9025 to 0, because G01 is the top-ranked "
            "gene. The port must read the GMT as tab-delimited.",
        mode=Mode.PARITY,
        gmt_description_with_spaces=True,
        nsets=8,
    ),
    Case(
        name="gmt_trailing_whitespace",
        why="`while(ss.good())` runs one iteration past the last token when a "
            "line ends in whitespace, and a failed `>>` into a std::string "
            "leaves it empty rather than untouched, so a phantom empty gene "
            "name joins the set. It never matches anything, but it inflates "
            "the Size column -- the number users read as the gene-set size -- "
            "by one. Confirmed against the binary. The port must not count it.",
        mode=Mode.PARITY,
        gmt_trailing_whitespace=True,
        nsets=8,
    ),
)

BY_NAME = {c.name: c for c in CASES}

PARITY_CASES = tuple(c for c in CASES if c.mode is Mode.PARITY)
CONTRACT_CASES = tuple(c for c in CASES if c.mode is Mode.CONTRACT)

#: Cases whose input file is deliberately malformed in a way that used to
#: change the helper's answer. Each was a ``Mode.DIVERGENCE`` case -- the port
#: was required to differ from a C++ that was defined and wrong -- until the
#: five bugs were fixed in the helper instead. They are ordinary parity cases
#: now, and the property they carry is stronger than before, because it is
#: required of *both* implementations: the answer on the malformed file must
#: equal the answer on its repair.
#:
#: Derived from the malformation flags rather than listed, so it cannot drift
#: from the cases themselves.
MALFORMED_CASES = tuple(
    c for c in CASES
    if (c.gmt_no_trailing_newline or c.rank_no_trailing_newline
        or c.gmt_blank_line or c.gmt_trailing_whitespace
        or c.gmt_description_with_spaces))

#: Cases a golden is recorded for: everything except the contract cases, where
#: the C++ reaches undefined behaviour and there is nothing to copy.
RECORDED_CASES = PARITY_CASES
