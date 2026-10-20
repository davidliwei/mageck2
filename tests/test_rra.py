"""Tests for the bundled RRA helper's own contract.

The smoke suite reaches RRA only through ``mageck2 test`` and ``mageck2
pathway``, and asserts little more than a zero exit status and a non-empty
summary file. Nothing pinned down RRA's own behaviour: how many permutations it
performs, what p value resolution that buys, or what that resolution costs in
achievable significance. This file covers that, so a change to the permutation
accounting has to state its effect rather than slip through.

The invariant worth protecting is the one in ``test_adjusted_p_floor_is_library
_size_independent``: the permutation count scales with the number of genes, and
``calculateFDRVector`` divides by that same number, so the smallest adjusted p
value a single strong hit can reach is fixed at about 0.005 whatever the library
size. A ceiling on the permutation count breaks that cancellation, which is the
cost ``--max-permutation-per-group`` trades run time for.

RRA is invoked as a subprocess by bare name, exactly as ``crisprFunction.py``
invokes it, so these tests exercise whichever build ``conftest.py`` put first on
PATH.
"""

import shutil
import subprocess

import numpy as np
import pytest

pytestmark = pytest.mark.skipif(
    shutil.which("RRA") is None,
    reason="the compiled RRA helper is not on PATH; run `pip install .` first",
)

#: Number of sgRNAs per gene in the fixtures. One value keeps every gene in a
#: single permutation group, which is what makes these runs quick enough for CI.
SGRNAS_PER_GENE = 4

#: Fraction of genes given a strongly depleting effect.
HIT_FRACTION = 0.02


def write_input(path, ngenes, seed=11):
    """Write an RRA input file: mostly null, with a planted set of hits.

    Format is the one ``fileio.cpp`` expects -- ``<item id> <group id>
    <list id> <value>`` behind a header row. Lower values rank first, so the
    planted genes are shifted down.
    """
    rng = np.random.default_rng(seed)
    hits = set(rng.choice(ngenes, size=max(1, int(ngenes * HIT_FRACTION)),
                          replace=False).tolist())

    rows = ["sgrna\tsymbol\tpool\tvalue"]
    for gene in range(ngenes):
        shift = -3.0 if gene in hits else 0.0
        for i in range(SGRNAS_PER_GENE):
            rows.append("sg%d_%d\tGENE%d\tlist\t%.6f"
                        % (gene, i, gene, rng.normal(shift, 1.0)))
    path.write_text("\n".join(rows) + "\n")
    return {"GENE%d" % g for g in hits}


def write_sparse_input(path, ngenes, nhits=1, seed=5):
    """Write an RRA input with a handful of unmistakable hits and nothing else.

    The fixtures above plant 2% hits, which hides what a permutation ceiling
    does to the achievable FDR: hundreds of genes tie at the p value floor, and
    the step-down pass in ``calculateFDRVector`` pulls the top gene's adjusted p
    value down by the size of that tied block. With a single hit there is no
    block, and the ``p * ngenes`` factor stands on its own. The shift is large
    enough (-6 sigma over every sgRNA) that the hit beats the whole null
    distribution, so its p value is exactly the floor ``1/(2N)``.
    """
    rng = np.random.default_rng(seed)
    rows = ["sgrna\tsymbol\tpool\tvalue"]
    for gene in range(ngenes):
        hit = gene < nhits
        name = "HIT%d" % gene if hit else "GENE%d" % gene
        for i in range(SGRNAS_PER_GENE):
            rows.append("sg%d_%d\t%s\tlist\t%.6f"
                        % (gene, i, name, rng.normal(-6.0 if hit else 0.0, 1.0)))
    path.write_text("\n".join(rows) + "\n")
    return {"HIT%d" % g for g in range(nhits)}


def run_rra(tmp_path, infile, *extra):
    """Run RRA and return (stdout, {gene: (lo_value, p, fdr)})."""
    outfile = tmp_path / "rra_out.txt"
    result = subprocess.run(
        ["RRA", "-i", str(infile), "-o", str(outfile),
         "-p", "0.25", "--skip-gene", "NA", *extra],
        cwd=tmp_path, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr

    rows = {}
    for line in outfile.read_text().splitlines()[1:]:
        f = line.split("\t")
        rows[f[0]] = (float(f[2]), float(f[3]), float(f[4]))
    assert rows, "RRA produced no gene rows"
    return result.stdout, rows


def announced_permutations(stdout):
    """The permutation count RRA reports capping to, or None if it did not."""
    for line in stdout.splitlines():
        if line.startswith("Capping permutations at "):
            return int(line.split()[3])
    return None


def implied_permutations(rows):
    """Recover the permutation count from the finest p value RRA could report.

    ``ComputePermutationbyNumSG`` assigns p as ``(i+j+1)/2/N`` over the sorted
    null, so a gene that beats the whole null distribution gets exactly
    ``1/(2N)``. That makes the permutation count observable from the output
    alone, which is what lets the tests below check it without a debug build.

    Returned as a float and compared approximately: ``fileio.cpp:296`` prints p
    with ``%10.4e``, so N is only recoverable to five significant digits.
    """
    return 1.0 / (2.0 * min(p for _, p, _ in rows.values()))


def test_rra_recovers_planted_hits(tmp_path):
    """The basic contract: strongly depleted genes rank first and clear FDR.

    Not specific to the permutation count -- this is the assertion that gives
    the rest of the file meaning, since a cap that broke the ranking would
    otherwise still satisfy every count-based check here.
    """
    infile = tmp_path / "in.txt"
    hits = write_input(infile, 2000)

    _, rows = run_rra(tmp_path, infile)

    ranked = sorted(rows, key=lambda g: rows[g][0])
    top = set(ranked[: len(hits)])
    assert len(top & hits) / len(hits) >= 0.9, (
        "planted hits did not dominate the top of the ranking: "
        "%d of %d recovered" % (len(top & hits), len(hits))
    )

    # A fraction, not every gene: four draws from N(-3,1) occasionally land
    # weakly, so demanding that every planted hit clears a threshold would be
    # asserting the fixture's luck. Recall at FDR 0.25 runs 0.97-1.00 across
    # seeds, so 0.9 has room without being vacuous.
    recall = sum(1 for g in hits if rows[g][2] <= 0.25) / len(hits)
    assert recall >= 0.9, "only %.2f of planted hits reached FDR 0.25" % recall


def test_p_values_resolve_to_the_permutation_count(tmp_path):
    """The finest reported p value is 1/(2N) for N permutations."""
    infile = tmp_path / "in.txt"
    write_input(infile, 2000)

    _, rows = run_rra(tmp_path, infile, "--permutation", "10")

    # rand_passnum=10 -> ComputePermutationbyNumSG runs (10+1)*2000 of them.
    assert implied_permutations(rows) == pytest.approx(11 * 2000, rel=1e-3)


def test_small_library_raises_the_permutation_count(tmp_path):
    """RRA.cpp:243 floors the count so small libraries still resolve p."""
    infile = tmp_path / "in.txt"
    write_input(infile, 200)

    stdout, rows = run_rra(tmp_path, infile)

    assert "Increase the number of permutations" in stdout
    assert announced_permutations(stdout) is None, "the floor must not also cap"
    assert implied_permutations(rows) >= 100000 * (1 - 1e-3)


def test_moderate_library_is_left_alone(tmp_path):
    """Between the floor and the ceiling, neither adjustment fires."""
    infile = tmp_path / "in.txt"
    write_input(infile, 2000)

    stdout, rows = run_rra(tmp_path, infile)

    assert "Increase the number of permutations" not in stdout
    assert announced_permutations(stdout) is None
    assert implied_permutations(rows) == pytest.approx(101 * 2000, rel=1e-3)


def test_a_large_library_is_not_capped_by_default(tmp_path):
    """Nothing bounds the count unless the caller asks for a bound.

    An 8,000-gene library permutes 101*8,000 = 808,000 times per sgRNA-count
    group and a genome-wide one 2.02 million. That looks like more resolution
    than an FDR threshold could use, and it would be if adjusted p values were
    the raw ones; the test below shows why it is not.
    """
    infile = tmp_path / "in.txt"
    write_input(infile, 8000)

    stdout, rows = run_rra(tmp_path, infile)

    assert announced_permutations(stdout) is None, (
        "the permutation count was bounded without --max-permutation-per-group"
    )
    assert implied_permutations(rows) == pytest.approx(101 * 8000, rel=1e-3)


@pytest.mark.parametrize("ngenes", [2000, 8000])
def test_adjusted_p_floor_is_library_size_independent(tmp_path, ngenes):
    """The invariant the permutation count exists to hold.

    ``calculateFDRVector`` multiplies the smallest p value by the number of
    genes, and the permutation count is ``(rand_passnum+1) * groupNum``, so for
    a gene that beats the entire null distribution the adjusted p value is
    ``groupNum / (2 * (rand_passnum+1) * groupNum)``: one over twice the pass
    count, with the library size cancelled out. At the default rand_passnum of
    100 that is 1/202, about 0.005, and a fourfold change in library size must
    not move it. Reducing the permutation count without reducing it in
    proportion to the gene count breaks this, and the break is invisible in
    every other test in this file.
    """
    infile = tmp_path / "in.txt"
    hits = write_sparse_input(infile, ngenes)

    _, rows = run_rra(tmp_path, infile)

    hit = hits.pop()
    assert rows[hit][2] == pytest.approx(1.0 / 202.0, rel=1e-2), (
        "the single planted hit reached FDR %.6f, not the expected 1/202"
        % rows[hit][2]
    )
    assert rows[hit][2] < 0.01, "an isolated hit must stay callable at FDR 0.01"


def test_a_ceiling_raises_the_adjusted_p_floor(tmp_path):
    """What a ceiling costs, stated rather than discovered later.

    Fixing the permutation count breaks the cancellation above, so the floor
    becomes ``ngenes / (2 * ceiling)`` and grows with the library. Asserted here
    so that the option's documented tradeoff is the measured one: 4,000 genes
    against a 20,000 ceiling puts the floor at 0.1, and the hit that reaches
    0.005 by default cannot be called at FDR 0.05 at all.
    """
    infile = tmp_path / "in.txt"
    hits = write_sparse_input(infile, 4000)

    _, rows = run_rra(tmp_path, infile, "--max-permutation-per-group", "20000")

    hit = hits.pop()
    assert rows[hit][2] == pytest.approx(4000 / (2.0 * 20000), rel=1e-2)
    assert rows[hit][2] > 0.05, (
        "this ceiling is supposed to cost the hit its FDR 0.05 call; it got "
        "%.4f, so the arithmetic above no longer describes the code"
        % rows[hit][2]
    )


@pytest.mark.parametrize("ngenes,ceiling", [(2000, 200000), (5000, 500000)])
def test_the_ceiling_is_never_exceeded(tmp_path, ngenes, ceiling):
    """The count compared against the ceiling must be the one actually run.

    The kernels run ``(rand_passnum+1) * groupNum`` permutations, so a guard
    written against ``rand_passnum * groupNum`` lets exactly one extra pass
    through: these two cases ran 202,004 and 504,999 against ceilings of
    200,000 and 500,000, and the announced number disagreed with the p values.
    """
    infile = tmp_path / "in.txt"
    write_input(infile, ngenes)

    stdout, rows = run_rra(tmp_path, infile,
                           "--max-permutation-per-group", str(ceiling))

    announced = announced_permutations(stdout)
    assert announced is not None, "the ceiling did not fire at its own boundary"
    assert announced <= ceiling, (
        "announced %d permutations against a ceiling of %d" % (announced, ceiling)
    )
    assert implied_permutations(rows) == pytest.approx(announced, rel=1e-3), (
        "RRA announced %d permutations but the p values imply %.0f"
        % (announced, implied_permutations(rows))
    )


def test_explicit_permutation_overrides_the_cap(tmp_path):
    """An explicit --permutation is honoured, however large it makes the count."""
    infile = tmp_path / "in.txt"
    write_input(infile, 8000)

    stdout, rows = run_rra(tmp_path, infile, "--permutation", "100",
                           "--max-permutation-per-group", "20000")

    assert announced_permutations(stdout) is None, (
        "--permutation was given, so RRA must not override it"
    )
    assert implied_permutations(rows) == pytest.approx(101 * 8000, rel=1e-3)


def test_max_permutation_per_group_is_honoured(tmp_path):
    """The ceiling is adjustable, and lowering it lowers the p resolution."""
    infile = tmp_path / "in.txt"
    write_input(infile, 2000)

    stdout, rows = run_rra(tmp_path, infile,
                           "--max-permutation-per-group", "20000")

    assert announced_permutations(stdout) == 20000
    assert implied_permutations(rows) == pytest.approx(20000, rel=1e-3)


def test_ceiling_below_the_gene_count_reports_what_it_actually_did(tmp_path):
    """A ceiling the kernels cannot express is raised, and said out loud.

    Permutations run in whole multiples of the gene count, so 2*groupNum is the
    floor. The announced number must be the count used rather than the one
    asked for, or the p values would silently disagree with the message.
    """
    infile = tmp_path / "in.txt"
    write_input(infile, 2000)

    stdout, rows = run_rra(tmp_path, infile, "--max-permutation-per-group", "100")

    assert announced_permutations(stdout) == 2 * 2000
    assert implied_permutations(rows) == pytest.approx(2 * 2000, rel=1e-3)


def test_capping_does_not_change_the_called_gene_set_when_hits_are_dense(tmp_path):
    """A ceiling leaves a hit-rich screen's gene set alone -- and only that one.

    This is the case a ceiling is defensible in, and the reason it is offered as
    an option. With 2% of genes planted as hits, a large block ties at the p
    value floor and the step-down pass keeps the adjusted values well clear of
    it, so a fourfold cut in permutations barely moves the called set. The
    sparse tests above are what stops this result from being read as a general
    one.

    Permutation p values carry Monte Carlo noise, so this compares gene sets
    rather than values, and the yardstick is that noise: two runs of the *same*
    settings under different seeds agree to about 0.98 at FDR 0.01, so a
    threshold much above that would be testing the seed. The fourfold cut here
    measures 0.986; eightfold measures 0.949, i.e. it leaves the noise floor,
    which is why the ceiling is the caller's choice and not a default.
    """
    infile = tmp_path / "in.txt"
    write_input(infile, 8000)

    _, capped = run_rra(tmp_path, infile,
                        "--max-permutation-per-group", "200000")
    _, uncapped = run_rra(tmp_path, infile)

    for threshold, min_jaccard in ((0.01, 0.95), (0.05, 0.95), (0.25, 0.95)):
        a = {g for g, (_, _, fdr) in capped.items() if fdr <= threshold}
        b = {g for g, (_, _, fdr) in uncapped.items() if fdr <= threshold}
        assert a and b, "no genes called at FDR %s" % threshold
        jaccard = len(a & b) / len(a | b)
        assert jaccard >= min_jaccard, (
            "capping changed the FDR<=%s gene set too much: Jaccard %.3f "
            "(%d capped, %d uncapped)" % (threshold, jaccard, len(a), len(b))
        )
