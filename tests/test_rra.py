"""Tests for the bundled RRA helper's own contract.

The smoke suite reaches RRA only through ``mageck2 test`` and ``mageck2
pathway``, and asserts little more than a zero exit status and a non-empty
summary file. Nothing pinned down RRA's own behaviour: how many permutations it
performs, what p value resolution that buys, or whether the gene set it calls
survives a change to that count. This file covers that, so a change to the
permutation accounting has to state its effect rather than slip through.

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


def test_large_library_caps_the_permutation_count(tmp_path):
    """Above the ceiling the count stops tracking the number of genes.

    Without the cap this library would permute 101*8000 = 808,000 times per
    sgRNA-count group, and a genome-wide one 2.02 million, resolving p far
    past anything an FDR threshold can use.
    """
    infile = tmp_path / "in.txt"
    write_input(infile, 8000)

    stdout, rows = run_rra(tmp_path, infile)

    announced = announced_permutations(stdout)
    assert announced is not None, "the cap did not fire on an 8,000-gene library"
    assert announced <= 500000
    assert announced < 101 * 8000, "capping left the count unchanged"
    assert implied_permutations(rows) == pytest.approx(announced, rel=1e-3), (
        "RRA announced %d permutations but the p values imply %.0f"
        % (announced, implied_permutations(rows))
    )


def test_explicit_permutation_overrides_the_cap(tmp_path):
    """An explicit --permutation is honoured, however large it makes the count."""
    infile = tmp_path / "in.txt"
    write_input(infile, 8000)

    stdout, rows = run_rra(tmp_path, infile, "--permutation", "100")

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


def test_capping_does_not_change_the_called_gene_set(tmp_path):
    """The point of the ceiling: fewer permutations, same genes.

    Permutation p values carry Monte Carlo noise, so this compares gene sets
    rather than values. The threshold is deliberately loose at FDR 0.01, where
    two runs of the *same* settings under different seeds agree to about 0.98 --
    tightening it would be testing the seed, not the cap.
    """
    infile = tmp_path / "in.txt"
    write_input(infile, 8000)

    _, capped = run_rra(tmp_path, infile)
    _, uncapped = run_rra(tmp_path, infile, "--permutation", "100")

    for threshold, min_jaccard in ((0.01, 0.95), (0.05, 0.95), (0.25, 0.95)):
        a = {g for g, (_, _, fdr) in capped.items() if fdr <= threshold}
        b = {g for g, (_, _, fdr) in uncapped.items() if fdr <= threshold}
        assert a and b, "no genes called at FDR %s" % threshold
        jaccard = len(a & b) / len(a | b)
        assert jaccard >= min_jaccard, (
            "capping changed the FDR<=%s gene set too much: Jaccard %.3f "
            "(%d capped, %d uncapped)" % (threshold, jaccard, len(a), len(b))
        )
