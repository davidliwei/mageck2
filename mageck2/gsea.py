"""Gene Set Enrichment Analysis, in Python.

This is the port of the bundled ``mageckGSEA`` C++ helper
(``gsea/src/calculategsea.cpp``). It exists for the same reason ``rra.py``
does: so that pathway enrichment can be called in process -- a mapping of
gene sets and two parallel sequences in, an object out -- instead of through
``subprocess`` with a temporary file on each side, which is what
``pathwayFunc.py`` does today.

The port is held to agreement with the C++ rather than to its own internal
consistency. ``tests/parity/gsea/`` records goldens from the compiled helper
and ``tests/test_gsea_parity.py`` compares this module against them column by
column. Several things here are therefore deliberately literal:

* The enrichment score is truncated to 32 bits before the permutation test
  compares it (``gs_ret::gscore`` is a ``float`` while everything computing it
  is a ``double``). This is not cosmetic: measured across the corpus, the
  truncation changes a permutation decision about once in 7,000 permutations,
  and a port that skips it fails 12 of the 15 parity cases.
* Gene sets are visited in lexicographic order, because the C++ iterates a
  ``std::map<string, ...>`` and that order determines what the generator is
  asked for and when. Visiting them in any other order draws the same numbers
  for different sets.
* The permutation sample is rejection sampling, drawing one uniform per
  attempt and discarding the attempt when the index is already taken. The
  rejected draws are part of the sequence, so they cannot be skipped and the
  loop cannot be replaced by a vectorized sample without changing every
  subsequent permutation in the run.
* ``parserankfile`` divides by one more than the number of ranked genes: the
  ``getline`` that reports end of file has already incremented the counter.
  That looks like an off-by-one and is kept, because it is what the shipped
  helper does and what its users' published numbers were computed with. A
  mutation that divides by the gene count instead is rejected by 14 of the 15
  parity cases.

What the port deliberately does *not* reproduce are the five input-handling
bugs the harness classifies as ``Mode.DIVERGENCE``: a missing final newline
dropping the last line of either input file, a blank GMT line becoming a gene
set named with the empty string, whitespace tokenization letting a gene set's
free-text description contribute gene names, and trailing whitespace adding a
phantom member. Each is required to differ, and the requirement is expressed
as a property -- the answer on a malformed file must equal the answer on its
repair -- rather than as a recorded number.

The generator is imported from :mod:`mageck2.rra` rather than reimplemented.
Both helpers vendor the same ``rngs.cpp``, verified byte-identical, so there is
one Lehmer stream in this package and not two subtly different ones. Moving it
to a module of its own belongs with the removal of the C++, when neither port
has to stay diff-minimal against a published branch.
"""

import re

import numpy as np
from scipy.special import betainc

from mageck2.rra import MODULUS, LehmerStream

#: ``GSEA_DEFAULT_SEED`` in ``calculategsea.cpp``. ``PlantSeeds`` stores
#: ``x % MODULUS`` and ``mageckGSEA`` never calls ``SelectStream``, so only
#: stream 0 is ever used.
DEFAULT_SEED = 123456

#: The largest seed the helper accepts. It rejects anything outside
#: ``1..MODULUS-1`` rather than folding it, because ``PutSeed`` reduces modulo
#: MODULUS and a seed of 0 or MODULUS would leave the generator stuck at 0.
MAX_SEED = MODULUS - 1

#: ``--permutation``'s default in ``argsParser.py``.
DEFAULT_PERMUTATION = 1000

#: Columns of the helper's output file, in order.
HEADER = "Pathway\tSize\tES\tp\tp_permutation\tFDR\tRanking\tHits\tLFC"

#: ``gsmax``'s initial value in ``getgseascore_fast``. The maximum is taken
#: over at least one walk position whenever a set is not degenerate, and every
#: walk value is bounded below by ``-1``, so the sentinel never wins and
#: ``ndarray.max`` is equivalent to the C++ loop. Kept as a named constant so
#: that equivalence is checkable rather than assumed.
_GSMAX_SENTINEL = -1000000.0

#: How many uniforms to draw from the generator at a time while rejection
#: sampling. The stream is consumed strictly in order, so the block size is
#: invisible in the result; it only trades memory for Python-level calls.
_DRAW_BLOCK = 1 << 16

#: How many permutations to score at once. The enrichment score of one
#: permutation depends only on its own indices, so batching changes nothing
#: and bounds the temporary at ``_ES_BATCH * num`` elements rather than at
#: ``perms * num``.
_ES_BATCH = 4096

#: ``strtod``'s leading-numeric-prefix rule. ``atof`` yields 0.0 for a field
#: with no numeric prefix rather than raising, which is how a non-numeric
#: header survives ``--skip-header 0``.
_LEADING_FLOAT = re.compile(r"\s*[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")


def _lehmer(seed):
    """The generator, with the helper's own bounds check on the seed."""
    seed = int(seed)
    if not 1 <= seed <= MAX_SEED:
        raise ValueError(
            "seed must be between 1 and %d; got %r" % (MAX_SEED, seed))
    return LehmerStream(seed, block=_DRAW_BLOCK)


def _atof(text):
    """``atof``: parse a leading number, return 0.0 when there is none."""
    match = _LEADING_FLOAT.match(text)
    return float(match.group()) if match else 0.0


def _format(value):
    """``ostream << double`` with default flags: six significant digits.

    Verified against the helper's output on the whole corpus: ``"%.6g"``
    reproduces it exactly for these values, which is what lets a golden
    recorded from the C++ be compared to this port as text.
    """
    return "%.6g" % value


def _lines(text):
    """The lines of ``text``, with a final newline not implying a last line.

    ``str.split("\\n")`` leaves a trailing empty string for a file that ends
    in a newline and does not for one that does not, so a real unterminated
    last line survives here. The C++ loses it: the ``getline`` that consumes
    it also sets ``eofbit``, and both readers test for end of file before
    using what they just read. That is two of the five divergences, and this
    one line is the whole of the fix.
    """
    parts = text.split("\n")
    if parts and parts[-1] == "":
        parts.pop()
    return parts


def read_gmt(path):
    """Read a GMT file into ``{set name: [member, ...]}``.

    GMT is tab-delimited: set name, a free-text description, then members.
    The C++ tokenizes on whitespace instead, so the second and later words of
    a description are read as members -- usually matching nothing and merely
    inflating the reported ``Size``, but gene symbols do appear in prose
    descriptions, and one that matches is scored. Splitting on tabs is the
    fix, and it also removes the phantom empty member a line with trailing
    whitespace otherwise gains.

    Blank lines are skipped. The C++ stores one as a gene set named with the
    empty string, which then matches the default ``--pathway-name`` of ``""``,
    takes the single-set branch, and reports that one empty set while silently
    skipping every real one at exit 0. Nothing here can produce a set named
    ``""``, which is what makes that unreachable.

    A repeated set name keeps the last definition, as ``PM[name] = ...`` does.
    """
    with open(path) as handle:
        text = handle.read()

    sets = {}
    for raw in _lines(text):
        line = raw.rstrip("\r")
        if not line.strip():
            continue
        fields = line.split("\t")
        name = fields[0].strip()
        if not name:
            continue
        sets[name] = [f.strip() for f in fields[2:] if f.strip()]
    return sets


def read_rank(path, score_column=2, skip_header=True):
    """Read a ranked gene list into ``(names, values, ncount)``.

    ``score_column`` is 1-based over the columns after the gene name, matching
    the helper's ``-c``. ``0`` selects the line-number branch described below.

    ``ncount`` is returned separately because it is *not* ``len(names)``: the
    C++ increments its line counter once more on the ``getline`` that reports
    end of file, and the walk weights are divided by that larger number. It is
    returned rather than recomputed so there is one place where the quirk is
    written down.

    ``parserankfile`` has two unrelated scoring branches and ``-c`` chooses
    between them. With a column given, the genes are sorted by ``(value,
    name)`` ascending and each gene's walk weight comes from its *position*;
    the values themselves are then discarded. With ``-c 0`` there is no
    column: a gene's value is its 1-based line number, left at 0 for any line
    with fewer than two fields, and the weight comes from that. Production
    never uses ``-c 0`` -- ``pathwayFunc.py`` always passes a column -- and it
    had no test coverage at all before the parity corpus.
    """
    with open(path) as handle:
        text = handle.read()

    rows = _lines(text)
    if skip_header and rows:
        rows = rows[1:]

    provide_column = score_column > 0
    names = []
    values = []
    for lineno, raw in enumerate(rows, start=1):
        fields = raw.rstrip("\r").split("\t")
        if len(fields) == 1:
            # The helper tokenizes on whitespace, so a space-separated file
            # parses there and would not here. Fall back only when the tab
            # split found nothing to split on, which leaves genuinely
            # tab-delimited files -- including ones whose fields contain
            # spaces -- read as written.
            fields = raw.split()
        if not fields or not fields[0]:
            continue
        if provide_column and len(fields) < score_column + 1:
            raise ValueError(
                "%s: line %d has %d column(s), so --ranking-column %d has "
                "nothing to read" % (path, lineno, len(fields), score_column))
        names.append(fields[0])
        if provide_column:
            values.append(_atof(fields[score_column]))
        else:
            # The C++ sets this only for a line with two or more fields,
            # leaving 0 otherwise, and 0 is what an absent value means here.
            values.append(float(lineno) if len(fields) >= 2 else 0.0)

    if not names:
        raise ValueError("%s contains no ranked genes" % path)
    # The end-of-file getline the C++ counts anyway.
    return names, values, len(names) + 1


def _walk_weights(ncount, n):
    """``|i * 2 / ncount - 1|`` for ``i`` in ``0..n-1``.

    The weight of a gene depends only on where it sits in the ranked list, so
    this is the whole of the scoring once the list is ordered.
    """
    i = np.arange(n, dtype=np.float64)
    return np.abs(i * 2.0 / ncount - 1.0)


def _ranked_list(names, values, ncount, score_column, reverse):
    """Order the gene list and give each position its walk weight."""
    n = len(names)
    if score_column > 0:
        # sort by (value, name) ascending, exactly as the C++ sorts its
        # vector of pairs, with the name breaking value ties. Without the
        # name tie-break 3 of the 15 parity cases disagree.
        order = sorted(range(n), key=lambda i: (values[i], names[i]))
        ordered = [names[i] for i in order]
        weights = _walk_weights(ncount, n)
    else:
        ordered = list(names)
        weights = np.abs(np.asarray(values, dtype=np.float64) * 2.0 / ncount
                         - 1.0)

    if reverse:
        # ``-e`` reverses the names and not the weights (calculategsea.cpp
        # :237-239), so position i keeps position i's weight and acquires a
        # different gene. Not the symmetry the flag's name suggests, and
        # reproduced rather than corrected because it is well defined.
        ordered.reverse()
    return ordered, weights


def _gene_index(names):
    """``{name: position}``, last occurrence winning.

    ``nameindex[gname[i]] = i`` in a loop, so a gene listed twice resolves to
    its later position. ``dict`` assignment has the same rule.
    """
    index = {}
    for i, name in enumerate(names):
        index[name] = i
    return index


def _members_to_positions(members, gene_index):
    """Sorted, duplicate-free positions of ``members`` in the ranked list.

    Members absent from the ranked list are dropped. The C++ looks them up
    with ``count`` before ``[]``, so -- unlike the bug fixed in #46 -- an
    absent gene is not inserted at position 0 and treated as top-ranked.
    """
    seen = set()
    positions = []
    for name in members:
        position = gene_index.get(name)
        if position is not None and position not in seen:
            seen.add(position)
            positions.append(position)
    positions.sort()
    return positions


def _enrichment_batch(weights, index):
    """32-bit enrichment score and argmax for each row of ``index``.

    ``index`` is ``(rows, num)``, each row sorted ascending. Returns
    ``(scores, argmax)`` with ``scores`` already truncated to ``float32``,
    because that truncation is what the permutation test compares at.

    This is the vectorized form of ``getgseascore_fast``. The scalar loop
    accumulates ``sumincval`` and takes a running maximum of ``decvalue *
    (ni - i) + sumincval``; here the accumulation is a ``cumsum`` along the
    row and the maximum is taken over the row at once. ``argmax`` returns the
    first maximum, which is what ``if (gsvalue > gsmax)`` selects.

    A set whose weights sum to zero is degenerate and scores 0.0, as are a set
    larger than the ranked list and an empty one. Those are not edge cases to
    tidy away: ``whole_list_set`` and ``absent_set`` in the corpus exist
    because real GMT files contain both.
    """
    rows, num = index.shape
    n = weights.shape[0]
    scores = np.zeros(rows, dtype=np.float64)
    argmax = np.zeros(rows, dtype=np.int64)
    if num <= 0 or n <= num:
        return scores.astype(np.float32), argmax

    values = weights[index]
    # Both sums are taken in the C++'s order, left to right over the sorted
    # overlap, because floating-point addition is not associative and the
    # permutation test compares scores that can cancel to around 1e-16.
    # ``ndarray.sum`` is pairwise, which is more accurate and therefore wrong
    # here; ``cumsum`` is sequential, so its last column is the running total
    # the C++ computes.
    running = np.cumsum(values, axis=1)
    total = running[:, -1]
    live = total != 0.0
    if live.any():
        decrement = -1.0 / (n - float(num))
        position = np.arange(num, dtype=np.float64)
        # ``sumincval += gscore[ni] / sumcrv``: each term is divided before it
        # is accumulated. Accumulating first and dividing once is a different
        # computation, and the difference was measured to flip one permutation
        # decision per run on sets whose score cancels.
        increment = np.cumsum(values[live] / total[live][:, None], axis=1)
        # Left as a separate multiply and add. The C++ compiler contracts the
        # same expression into a fused multiply-add wherever FMA is in the
        # baseline ISA, which moves the score of a set whose weights cancel to
        # zero by about 5.6e-17 -- so the helper does not even agree with
        # itself across architectures there. The comparator allows ES that
        # much slack for exactly this reason, and the 32-bit truncation below
        # keeps the permutation counts exact regardless.
        walk = decrement * (index[live] - position) + increment
        scores[live] = walk.max(axis=1)
        argmax[live] = walk.argmax(axis=1)
    return scores.astype(np.float32), argmax


class _Draws:
    """The generator's output as one continuous sequence of gene positions.

    Blocks are an implementation detail of not calling into NumPy once per
    draw. The sequence they come from must stay continuous *across gene
    sets*, because the C++ never restarts it: the sample for one set begins
    wherever the previous one stopped. Keeping a block's unconsumed tail is
    therefore a correctness requirement and not an optimization -- discarding
    it on each set shifts every later set's permutations, which shows up as a
    handful of rank inversions in ``p_permutation`` and in nothing else at
    all, since the enrichment scores and their beta p values do not touch the
    generator.
    """

    def __init__(self, stream, nsz, block=_DRAW_BLOCK):
        self._stream = stream
        self._nsz = int(nsz)
        self._block = int(block)
        self._buf = []
        self._cursor = 0

    def _refill(self):
        draws = self._stream.uniform(self._block) * self._nsz
        # ``cid = (int)(Random() * nsz)`` then clamped to nsz - 1. The clamp
        # cannot fire -- Random() is strictly below 1 because the state stays
        # in 1..MODULUS-1 -- and is kept because the arithmetic is the C++'s.
        np.minimum(draws, self._nsz - 1, out=draws)
        self._buf = draws.astype(np.int64).tolist()
        self._cursor = 0

    def fill(self, taken, scratch):
        """Fill ``taken`` with rejection-sampled gene positions.

        ``taken`` is ``(perms, num)``; each row becomes a sample of ``num``
        distinct positions drawn without replacement, left sorted ascending.

        The draw order is the whole point. The C++ computes a candidate,
        clamps it, and discards it when that position is already in the
        sample -- so a rejected attempt still consumes a number. How many
        numbers a row consumes is not known until it is drawn, which is why
        this is a sequential walk and not a vectorized sample. Collisions are
        common rather than negligible: for 60 members out of 2,000 genes more
        than half of all rows hit one, so there is no fast path worth adding.

        ``scratch`` is a ``bytearray`` of length ``nsz`` marking taken
        positions, cleared per row so the cost does not scale with the gene
        list on every permutation.
        """
        perms, num = taken.shape
        buf = self._buf
        cursor = self._cursor
        size = len(buf)
        for row in range(perms):
            touched = []
            filled = 0
            while filled < num:
                if cursor >= size:
                    self._cursor = cursor
                    self._refill()
                    buf = self._buf
                    cursor = 0
                    size = len(buf)
                cid = buf[cursor]
                cursor += 1
                if not scratch[cid]:
                    scratch[cid] = 1
                    touched.append(cid)
                    taken[row, filled] = cid
                    filled += 1
            for cid in touched:
                scratch[cid] = 0
        self._cursor = cursor
        taken.sort(axis=1)
        return taken


def _permutation_p(draws, weights, num, perms, observed, scratch):
    """Fraction of permuted sets scoring above ``observed``.

    ``observed`` is compared as a 32-bit value against 32-bit permutation
    scores, which is what the C++ does by storing both in a ``float`` field.
    """
    if perms <= 0:
        return 0.0
    threshold = np.float64(np.float32(observed))
    above = 0
    done = 0
    while done < perms:
        rows = min(_ES_BATCH, perms - done)
        taken = np.empty((rows, num), dtype=np.int64)
        draws.fill(taken, scratch)
        scores, _ = _enrichment_batch(weights, taken)
        above += int(np.count_nonzero(scores.astype(np.float64) > threshold))
        done += rows
    return above * 1.0 / perms


class GSEAResult:
    """One row per gene set, in the order the helper writes them.

    ``to_text`` renders the historical file format, which is what lets the
    goldens recorded from the C++ go on serving as the oracle once the C++ is
    gone.
    """

    HEADER = HEADER

    def __init__(self, names, sizes, es, p, p_permutation, fdr, ranking,
                 hits):
        self.names = list(names)
        self.sizes = np.asarray(sizes, dtype=np.int64)
        self.es = np.asarray(es, dtype=np.float64)
        self.p = np.asarray(p, dtype=np.float64)
        self.p_permutation = np.asarray(p_permutation, dtype=np.float64)
        self.fdr = np.asarray(fdr, dtype=np.float64)
        #: 1-based position after sorting by p, or 0 for every row when the
        #: output is unsorted -- the C++ writes a literal 0 there
        #: (``calculategsea.cpp:494``). A column of identical values carries
        #: no ordering claim, and the comparator checks this column only when
        #: it holds more than one distinct value.
        self.ranking = np.asarray(ranking, dtype=np.int64)
        #: ``gindex``: the position *within the sorted overlap* at which the
        #: running enrichment peaks. Not a count of anything, despite the
        #: column's name.
        self.hits = np.asarray(hits, dtype=np.int64)

    def __len__(self):
        return len(self.names)

    def __iter__(self):
        for i, name in enumerate(self.names):
            yield {
                "pathway": name,
                "size": int(self.sizes[i]),
                "es": float(self.es[i]),
                "p": float(self.p[i]),
                "p_permutation": float(self.p_permutation[i]),
                "fdr": float(self.fdr[i]),
                "ranking": int(self.ranking[i]),
                "hits": int(self.hits[i]),
            }

    def to_text(self):
        """The helper's output format, tab-separated with a header row."""
        out = [self.HEADER]
        for i, name in enumerate(self.names):
            out.append("\t".join([
                name,
                str(int(self.sizes[i])),
                _format(float(self.es[i])),
                _format(float(self.p[i])),
                _format(float(self.p_permutation[i])),
                _format(float(self.fdr[i])),
                str(int(self.ranking[i])),
                str(int(self.hits[i])),
                # The helper writes a literal 0 in both output branches. It is
                # not a log fold change that happened to be zero.
                "0",
            ]))
        return "\n".join(out) + "\n"


def _benjamini_hochberg(p_values):
    """``calculateFDRVector``, transcribed.

    Ascending p values in, adjusted values out. Only the last entry is clamped
    to 1 -- the C++ clamps that one and then takes a running minimum backwards,
    which is what pulls the rest below 1 -- so clamping every entry instead
    would be a different function on inputs where it matters.
    """
    n = len(p_values)
    if n == 0:
        return []
    adjusted = [p_values[i] / (i + 1) * n for i in range(n)]
    if adjusted[-1] > 1.0:
        adjusted[-1] = 1.0
    for i in range(n - 2, -1, -1):
        if adjusted[i] > adjusted[i + 1]:
            adjusted[i] = adjusted[i + 1]
    return adjusted


def _score_one(name, members, ranked, weights, gene_index, draws, perms,
               scratch):
    """Score one gene set: ``(es, hits, p, p_permutation)``.

    Returns the degenerate answer -- score 0, p 1, permutation p 1 -- without
    touching the generator when the set has no usable overlap, which is what
    the C++ does by returning before its permutation loop. That matters beyond
    tidiness: a set that consumed draws here would shift every later set's
    permutations.
    """
    positions = _members_to_positions(members, gene_index)
    if not positions:
        return np.float32(0.0), 0, 1.0, 1.0

    index = np.asarray(positions, dtype=np.int64)[None, :]
    scores, argmax = _enrichment_batch(weights, index)
    es = scores[0]
    hits = int(argmax[0])
    num = index.shape[1]
    n = weights.shape[0]
    if n <= num or weights[index[0]].sum() == 0.0:
        return es, hits, 1.0, 1.0

    # BetaNoncentralCdf(a, b, 0.0, x, 1e-10) with lambda 0 at every call site,
    # which is scipy's regularized incomplete beta exactly.
    percentile = positions[hits] * 1.0 / float(n)
    p = float(betainc(hits + 1, num - hits, percentile))
    p_permutation = _permutation_p(draws, weights, num, perms, float(es),
                                   scratch)
    return es, hits, p, p_permutation


def enrich(pathways, genes, scores, perms=DEFAULT_PERMUTATION,
           seed=DEFAULT_SEED, sort_by_p=True, reverse=False,
           pathway_name=None, ncount=None):
    """Test each gene set in ``pathways`` for enrichment in a ranked list.

    ``pathways`` maps a set name to its members. ``genes`` and ``scores`` are
    parallel sequences: a gene and the statistic it is ranked on, smaller
    first. Only the *order* ``scores`` induces is used -- the walk weight of a
    gene comes from its position, not from its value -- so any monotone
    transform of ``scores`` gives the same answer.

    ``pathway_name`` restricts the run to one named set, which reports it
    without an FDR pass and with ``Ranking`` fixed at 1. A name that is not in
    ``pathways`` runs the full analysis, as the C++ does.

    ``ncount`` overrides the divisor of the walk weights. It exists so
    :func:`enrich_from_file` can pass the helper's line count, which is one
    more than the number of genes; leave it alone otherwise.

    Returns a :class:`GSEAResult`.
    """
    genes = [str(g) for g in genes]
    if len(genes) != len(scores):
        raise ValueError("genes and scores must be parallel")
    if not genes:
        raise ValueError("no ranked genes to test against")
    if ncount is None:
        ncount = len(genes) + 1

    ranked, weights = _ranked_list(genes, list(scores), ncount, 2, reverse)
    gene_index = _gene_index(ranked)
    # One sampler for the whole run: the generator is not restarted per set.
    draws = _Draws(_lehmer(seed), len(ranked))
    scratch = bytearray(len(ranked))

    if pathway_name is not None and pathway_name in pathways:
        members = pathways[pathway_name]
        es, hits, p, p_perm = _score_one(
            pathway_name, members, ranked, weights, gene_index, draws, perms,
            scratch)
        # The single-set branch reports the permutation p in the FDR column
        # and a literal 1 as the ranking; there is nothing to adjust against.
        return GSEAResult([pathway_name], [len(members)], [float(es)], [p],
                          [p_perm], [p_perm], [1], [hits])

    rows = []
    # Lexicographic, because the C++ walks a std::map and that fixes which
    # set gets which draws.
    for name in sorted(pathways):
        members = pathways[name]
        es, hits, p, p_perm = _score_one(
            name, members, ranked, weights, gene_index, draws, perms,
            scratch)
        rows.append([name, len(members), float(es), p, p_perm, hits])

    if not sort_by_p:
        return GSEAResult(
            [r[0] for r in rows], [r[1] for r in rows], [r[2] for r in rows],
            [r[3] for r in rows], [r[4] for r in rows],
            # Unsorted output repeats the permutation p in the FDR column and
            # writes 0 for every ranking.
            [r[4] for r in rows], [0] * len(rows), [r[5] for r in rows])

    # Ascending permutation p, with a larger enrichment score first on a tie.
    # The C++ sorts with an unstable quicksort, so sets tying on both keys may
    # come out in either order; the FDR pass is unaffected, because its
    # backward running minimum collapses a tied run to one value whatever the
    # order inside it.
    rows.sort(key=lambda r: (r[4], -r[2]))
    fdr = _benjamini_hochberg([r[4] for r in rows])
    return GSEAResult(
        [r[0] for r in rows], [r[1] for r in rows], [r[2] for r in rows],
        [r[3] for r in rows], [r[4] for r in rows], fdr,
        list(range(1, len(rows) + 1)), [r[5] for r in rows])


def enrich_from_file(gmt_path, rank_path, score_column=2,
                     perms=DEFAULT_PERMUTATION, seed=DEFAULT_SEED,
                     sort_by_p=True, skip_header=True, reverse=False,
                     pathway_name=None):
    """Run GSEA over a GMT file and a ranked gene list, mirroring the CLI.

    ``score_column`` is the helper's ``-c``: 1-based over the columns after
    the gene name, or 0 to rank by line number instead.
    """
    pathways = read_gmt(gmt_path)
    names, values, ncount = read_rank(
        rank_path, score_column=score_column, skip_header=skip_header)

    if score_column > 0:
        return enrich(pathways, names, values, perms=perms, seed=seed,
                      sort_by_p=sort_by_p, reverse=reverse,
                      pathway_name=pathway_name, ncount=ncount)

    # The line-number branch scores from the values themselves rather than
    # from the induced order, so it cannot go through enrich()'s sort.
    ranked, weights = _ranked_list(names, values, ncount, 0, reverse)
    gene_index = _gene_index(ranked)
    # One sampler for the whole run: the generator is not restarted per set.
    draws = _Draws(_lehmer(seed), len(ranked))
    scratch = bytearray(len(ranked))

    if pathway_name is not None and pathway_name in pathways:
        members = pathways[pathway_name]
        es, hits, p, p_perm = _score_one(
            pathway_name, members, ranked, weights, gene_index, draws, perms,
            scratch)
        return GSEAResult([pathway_name], [len(members)], [float(es)], [p],
                          [p_perm], [p_perm], [1], [hits])

    rows = []
    for name in sorted(pathways):
        members = pathways[name]
        es, hits, p, p_perm = _score_one(
            name, members, ranked, weights, gene_index, draws, perms,
            scratch)
        rows.append([name, len(members), float(es), p, p_perm, hits])

    if not sort_by_p:
        return GSEAResult(
            [r[0] for r in rows], [r[1] for r in rows], [r[2] for r in rows],
            [r[3] for r in rows], [r[4] for r in rows], [r[4] for r in rows],
            [0] * len(rows), [r[5] for r in rows])

    rows.sort(key=lambda r: (r[4], -r[2]))
    fdr = _benjamini_hochberg([r[4] for r in rows])
    return GSEAResult(
        [r[0] for r in rows], [r[1] for r in rows], [r[2] for r in rows],
        [r[3] for r in rows], [r[4] for r in rows], fdr,
        list(range(1, len(rows) + 1)), [r[5] for r in rows])
