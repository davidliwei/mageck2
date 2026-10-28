"""Tests for the bundled mageckGSEA helper's permutation contract.

``test_smoke.py`` reaches mageckGSEA through ``mageck2 pathway`` and through a
few direct invocations, and checks that scores and p-values are well formed.
Nothing pinned down whether a run is *repeatable*, and it was not:
``init_vec`` seeded the permutation test with ``srand(time(NULL))`` and drew
with ``rand()``, so two identical invocations more than a second apart returned
different permutation p-values, different FDRs, and a different pathway
ranking. Published pathway results could not be reproduced, including by the
person who produced them.

The permutations now run on the Lehmer generator in ``rngs.cpp`` -- the one RRA
already uses -- seeded from ``--seed``, which defaults to 123456. That makes a
run a function of its inputs, and makes Monte Carlo noise something a caller
can measure on purpose by varying the seed rather than something they get by
accident.

Only the helper is tested here. ``mageck2 pathway`` passes no seed and nothing
else on that path is random, so the helper's determinism is the command's.

mageckGSEA is invoked as a subprocess by bare name, exactly as
``pathwayFunc.py`` invokes it, so these tests exercise whichever build
``conftest.py`` put first on PATH.
"""

import shutil
import subprocess
import time

import numpy as np
import pytest

pytestmark = pytest.mark.skipif(
    shutil.which("mageckGSEA") is None,
    reason="the compiled mageckGSEA helper is not on PATH; run `pip install .` first",
)

#: The default ``--seed``, as declared in ``calculategsea.cpp``.
DEFAULT_SEED = 123456

#: One below MODULUS in ``rngs.cpp``; the largest seed the generator accepts.
MAX_SEED = 2147483646

#: Columns of the output that are computed from the input alone. None of them
#: may move when only the seed changes.
DETERMINISTIC_COLUMNS = ["Size", "ES", "p", "Hits", "LFC"]


def write_inputs(tmp_path, ngenes=400, nsets=12, setsize=25, seed=17):
    """Write a rank file and a GMT, sized so permutation p-values are informative.

    The sets have to be small against the ranked list -- ``getgseascore_fast``
    treats a set covering the whole list as degenerate and returns p=1 without
    permuting, which would make every assertion below vacuous.
    """
    rng = np.random.default_rng(seed)
    genes = ["GENE%d" % i for i in range(ngenes)]

    rank = tmp_path / "rank.txt"
    rank.write_text(
        "id\tscore\n"
        + "".join("%s\t%.6f\n" % (g, s) for g, s in zip(genes, rng.random(ngenes)))
    )

    gmt = tmp_path / "sets.gmt"
    gmt.write_text(
        "".join(
            "PATH%d\tna\t%s\n"
            % (i, "\t".join(rng.choice(genes, size=setsize, replace=False)))
            for i in range(nsets)
        )
    )
    return rank, gmt


def run_gsea(tmp_path, rank, gmt, *extra, tag="out", perms=200, check=True):
    """Run mageckGSEA and return (CompletedProcess, raw output text).

    ``tag`` only names the output file, so two runs in one test can be
    compared without one overwriting the other.
    """
    out = tmp_path / ("%s.txt" % tag)
    result = subprocess.run(
        ["mageckGSEA", "-H", "-s", "-c", "1", "-p", str(perms),
         "-g", str(gmt), "-r", str(rank), "-o", str(out), *extra],
        cwd=tmp_path, capture_output=True, text=True,
    )
    if check:
        assert result.returncode == 0, result.stderr
        return result, out.read_text()
    return result, None


def parse(text):
    """Parse mageckGSEA output into {pathway: {column: value}}."""
    lines = text.splitlines()
    header = lines[0].split("\t")
    rows = {}
    for line in lines[1:]:
        f = line.split("\t")
        rows[f[0]] = dict(zip(header[1:], f[1:]))
    assert rows, "mageckGSEA produced no pathway rows"
    return rows


def test_identical_runs_produce_identical_output(tmp_path):
    """The regression test: a run must be a function of its inputs.

    The sleep is the point of the test, not incidental. The old seeding read
    the clock in whole seconds, so two runs inside the same second agreed even
    when the bug was present -- a back-to-back comparison would have passed
    against ``srand(time(NULL))`` and caught nothing. Letting the clock tick
    over is what makes this a regression test.
    """
    rank, gmt = write_inputs(tmp_path)

    _, first = run_gsea(tmp_path, rank, gmt, tag="first")
    time.sleep(1.1)
    _, second = run_gsea(tmp_path, rank, gmt, tag="second")

    assert first == second, (
        "two identical runs disagreed, so the permutation test is still "
        "seeded from something outside its inputs"
    )


def test_the_default_seed_is_the_documented_one(tmp_path):
    """Passing the documented default explicitly must change nothing."""
    rank, gmt = write_inputs(tmp_path)

    _, implicit = run_gsea(tmp_path, rank, gmt, tag="implicit")
    _, explicit = run_gsea(tmp_path, rank, gmt, "--seed", str(DEFAULT_SEED),
                           tag="explicit")

    assert implicit == explicit


def test_a_different_seed_moves_only_the_monte_carlo_columns(tmp_path):
    """A seed change must move the permutation p-values and nothing else.

    This is what makes the seed useful rather than merely present: the
    enrichment score and the analytic p-value are properties of the input, so
    varying the seed measures the permutation test's own noise without
    disturbing them. It is also the check that would fail if the seed were
    accidentally wired into the scoring path.
    """
    rank, gmt = write_inputs(tmp_path)

    _, a = run_gsea(tmp_path, rank, gmt, tag="default_seed")
    _, b = run_gsea(tmp_path, rank, gmt, "--seed", "999", tag="other_seed")
    A, B = parse(a), parse(b)

    assert set(A) == set(B)
    for column in DETERMINISTIC_COLUMNS:
        differing = [p for p in A if A[p][column] != B[p][column]]
        assert not differing, (
            "%s depends on the permutation seed, which it must not: %s"
            % (column, differing[:5])
        )

    moved = [p for p in A if A[p]["p_permutation"] != B[p]["p_permutation"]]
    assert moved, (
        "no permutation p-value moved when the seed changed, so --seed is "
        "not reaching the permutation test"
    )


@pytest.mark.parametrize("seed", [0, -1, MAX_SEED + 1])
def test_an_out_of_range_seed_is_rejected(tmp_path, seed):
    """Seeds the generator cannot use must be refused, not silently accepted.

    ``rngs.cpp`` keeps its state in 1..MODULUS-1. ``PutSeed`` reduces a
    positive seed modulo MODULUS, so MODULUS itself becomes 0 -- a fixed point
    that makes ``Random()`` return 0 forever and every permutation draw pick
    gene 0. A negative seed sends ``PutSeed`` to the clock, and a seed of 0
    sends it to an interactive ``scanf`` prompt, which would hang a helper that
    is only ever run as a subprocess.
    """
    rank, gmt = write_inputs(tmp_path)

    result, _ = run_gsea(tmp_path, rank, gmt, "--seed", str(seed),
                         perms=10, check=False)

    assert result.returncode != 0, "an unusable seed was accepted"
    # The message, not just the exit status: a build without the option at all
    # also exits non-zero, with TCLAP's "couldn't find match for argument", so
    # checking only the status would let this test pass against the very code
    # it exists to rule out.
    assert "must be between" in result.stderr, (
        "the failure did not report the accepted seed range: %r" % result.stderr
    )


@pytest.mark.parametrize("seed", [1, MAX_SEED])
def test_the_seed_range_boundaries_are_accepted(tmp_path, seed):
    """The ends of the documented range must work, or the range is wrong."""
    rank, gmt = write_inputs(tmp_path)

    result, text = run_gsea(tmp_path, rank, gmt, "--seed", str(seed), perms=10)

    assert result.returncode == 0
    assert parse(text)


def test_the_seed_is_reported(tmp_path):
    """The run must say which seed it used, so a result can be reproduced.

    A default that is never printed is only half a fix: anyone holding an old
    output file needs the log to tell them which stream produced it.
    """
    rank, gmt = write_inputs(tmp_path)

    result, _ = run_gsea(tmp_path, rank, gmt, "--seed", "4242", perms=10)

    assert "4242" in result.stderr, (
        "the seed was not reported in the run log: %r" % result.stderr
    )
