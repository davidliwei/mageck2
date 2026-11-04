"""Tests for the Python RRA implementation and the switch that selects it.

Agreement with the compiled helper is established in ``test_rra_parity.py``
against recorded goldens. This file covers what that cannot: the public array
API, the behaviours the port *chooses* rather than inherits, the correctness of
the one optimization in it, and the environment switch itself.

The switch is tested because a migration switch that silently does nothing is
the worst possible outcome -- every parity test would still pass, against the
old implementation, and the port would appear validated without ever having
run.
"""

import os
import shutil
import subprocess
import sys

import numpy as np
import pytest

from mageck2 import crisprFunction, rra


# --------------------------------------------------------------------------
# The migration switch
# --------------------------------------------------------------------------

def test_the_python_implementation_is_the_default(monkeypatch):
    """The port is now what runs unless the environment says otherwise.

    This assertion is the switch's whole point, so it is worth being explicit
    about what it protects: if the default silently reverted to the compiled
    helper, every parity test would still pass -- against the old
    implementation -- and nothing else in the suite would notice.
    """
    monkeypatch.delenv(crisprFunction.RRA_IMPL_ENV, raising=False)
    assert crisprFunction.rra_implementation() == "python"


@pytest.mark.parametrize("value", ["cpp", "CPP", " cpp "])
def test_the_compiled_helper_remains_reachable(monkeypatch, value):
    """The fallback has to keep working while the C++ is still in the tree.

    It is what makes flipping the default reversible without a release.
    """
    monkeypatch.setenv(crisprFunction.RRA_IMPL_ENV, value)
    assert crisprFunction.rra_implementation() == "cpp"


@pytest.mark.parametrize("value", ["python", "Python", " PYTHON "])
def test_the_switch_selects_the_port(monkeypatch, value):
    monkeypatch.setenv(crisprFunction.RRA_IMPL_ENV, value)
    assert crisprFunction.rra_implementation() == "python"


def test_an_unknown_implementation_is_refused(monkeypatch):
    """Rather than falling back to a default and running the wrong code."""
    monkeypatch.setenv(crisprFunction.RRA_IMPL_ENV, "rust")
    with pytest.raises(SystemExit):
        crisprFunction.rra_implementation()


# --------------------------------------------------------------------------
# The random number generator
# --------------------------------------------------------------------------

def test_the_generator_matches_a_scalar_reference_across_block_boundaries():
    """The blockwise jump must agree with the plain recurrence.

    ``LehmerStream`` draws a block at a time by jumping ahead rather than
    iterating, which is what keeps the permutation null from being slower than
    the C++ it replaces. That makes the block boundary the one place it could
    be subtly wrong -- a reused or skipped state would be invisible in
    aggregate and fatal to the parity comparison -- so the block size here is
    deliberately small and coprime to the request size.
    """
    seed = rra.DEFAULT_SEED
    state = seed
    expected = []
    for _ in range(50):
        state = (rra.MULTIPLIER * state) % rra.MODULUS
        expected.append(state / rra.MODULUS)

    got = rra.LehmerStream(seed, block=7).uniform(50)
    assert np.array_equal(got, np.asarray(expected))


def test_the_generator_is_continuous_across_separate_calls():
    """Two calls must not restart or skip the stream."""
    one = rra.LehmerStream(rra.DEFAULT_SEED, block=7)
    first = np.concatenate([one.uniform(13), one.uniform(37)])
    whole = rra.LehmerStream(rra.DEFAULT_SEED, block=7).uniform(50)
    assert np.array_equal(first, whole)


def test_the_generator_stays_inside_the_open_unit_interval():
    """Nothing may return exactly 0 or 1.

    The control resampling multiplies by the number of control sgRNAs and
    truncates, so a draw of exactly 1.0 would index one past the end. The C++
    guards that with a clamp; the guarantee that makes the clamp unreachable is
    this property, and it is worth asserting rather than assuming.
    """
    values = rra.LehmerStream(rra.DEFAULT_SEED).uniform(200000)
    assert values.min() > 0.0
    assert values.max() < 1.0


@pytest.mark.parametrize("seed", [0, -1, rra.MODULUS])
def test_an_out_of_range_seed_is_refused(seed):
    with pytest.raises(ValueError):
        rra.LehmerStream(seed)


# --------------------------------------------------------------------------
# The nearest-element search
# --------------------------------------------------------------------------

def test_the_search_returns_a_present_element_not_an_insertion_point():
    """This is the distinction that makes ``searchsorted`` the wrong function.

    ``bTreeSearchingF`` answers "which element is closest", so every result
    indexes a value that exists. ``searchsorted`` answers "where would this go",
    which is a different question with a systematically different answer.
    """
    a = np.array([0.0, 1.0, 2.0, 3.0, 10.0])
    assert rra.nearest_index(a, [2.1])[0] == 2
    assert rra.nearest_index(a, [2.9])[0] == 3
    assert rra.nearest_index(a, [-5.0])[0] == 0
    assert rra.nearest_index(a, [99.0])[0] == 4
    # Nearest, not bracketing: 5.0 lies between 3 and 10 but is closer to 3.
    assert rra.nearest_index(a, [5.0])[0] == 3
    # Exactly equidistant between two neighbours: the C++ resolves the
    # two-element window with `fabs(a[hi]-v) > fabs(a[lo]-v)`, which is false
    # on a tie, so the upper index wins.
    assert rra.nearest_index(a, [2.5])[0] == 3


def test_the_search_is_vectorized_consistently():
    rng = np.random.default_rng(3)
    a = np.sort(rng.normal(size=500))
    q = rng.normal(size=200)
    batch = rra.nearest_index(a, q)
    one_at_a_time = [rra.nearest_index(a, [x])[0] for x in q]
    assert np.array_equal(batch, np.asarray(one_at_a_time))


# --------------------------------------------------------------------------
# The array API
# --------------------------------------------------------------------------

def _tiny_screen(ngenes=40, per_gene=4, seed=5):
    rng = np.random.default_rng(seed)
    items, genes, values = [], [], []
    for g in range(ngenes):
        shift = -4.0 if g < 2 else 0.0
        for i in range(per_gene):
            items.append("sg%d_%d" % (g, i))
            genes.append("GENE%02d" % g)
            values.append(rng.normal(shift, 1.0))
    return items, genes, np.asarray(values)


def test_the_api_takes_arrays_and_returns_an_object():
    """The point of the port: no file path in, no file path out."""
    items, genes, values = _tiny_screen()
    result = rra.rank(items, genes, values, alpha=0.25)

    assert len(result) == 40
    assert result.lo_value.shape == (40,)
    assert result.p.shape == (40,)
    assert np.all(result.p > 0) and np.all(result.p <= 1)
    assert np.all(result.fdr <= 1.0)
    assert result.goodsgrna.dtype == np.int64
    # The planted genes are the two strongest.
    ranked = [result.names[i] for i in np.argsort(result.lo_value)]
    assert set(ranked[:2]) == {"GENE00", "GENE01"}


def test_the_api_is_reproducible_without_being_asked():
    items, genes, values = _tiny_screen()
    a = rra.rank(items, genes, values).to_text()
    b = rra.rank(items, genes, values).to_text()
    assert a == b


def test_a_different_seed_moves_the_p_values_but_not_the_lo_values():
    """lo_value is computed from the data; p from the permutation null."""
    items, genes, values = _tiny_screen()
    a = rra.rank(items, genes, values, seed=123456)
    b = rra.rank(items, genes, values, seed=99999)
    assert np.array_equal(a.lo_value, b.lo_value)
    assert not np.array_equal(a.p, b.p)


def test_mismatched_input_lengths_are_refused():
    items, genes, values = _tiny_screen()
    with pytest.raises(ValueError):
        rra.rank(items[:-1], genes, values)


def test_a_control_list_matching_nothing_is_refused():
    """The contract chosen for issue #51's out-of-bounds read.

    Falling back to a uniform null would report p values that are not the
    control-based ones the caller asked for, and nothing in the output would
    say so. The message has to name the control list, because that is the input
    the caller has to fix.
    """
    items, genes, values = _tiny_screen()
    with pytest.raises(ValueError) as excinfo:
        rra.rank(items, genes, values, control=["not_in_this_library"])
    assert "control" in str(excinfo.value).lower()


def test_a_control_list_is_usable_when_it_matches():
    items, genes, values = _tiny_screen()
    controls = [s for s in items if s.startswith(("sg30", "sg31", "sg32"))]
    assert controls
    result = rra.rank(items, genes, values, control=controls)
    # Genes made up entirely of control sgRNAs are suppressed from the output,
    # as fileio.cpp:331 does.
    assert result.suppressed.sum() == 3
    assert "GENE30" not in result.to_text()


def test_blank_lines_in_a_control_file_are_ignored(tmp_path):
    """The C++ reader turns the trailing newline into an empty-named entry.

    It inflates the count it reports and makes it warn about an sgRNA with no
    name. Nothing downstream is affected -- an empty name matches no sgRNA --
    so the port drops blanks instead of reproducing the artifact.
    """
    path = tmp_path / "control.txt"
    path.write_text("sg1\n\nsg2\n   \n")
    assert rra.read_control(str(path)) == ["sg1", "sg2"]


# --------------------------------------------------------------------------
# Both implementations through the real CLI
# --------------------------------------------------------------------------

@pytest.mark.skipif(shutil.which("mageck2") is None or shutil.which("RRA") is None,
                    reason="needs an installed mageck2 and the compiled helper")
def test_both_implementations_give_the_same_gene_summary(tmp_path):
    """The end-to-end check: a real `mageck2 test` run, both ways.

    The parity suite compares RRA's own output file. This compares what a user
    actually reads, after normalization, LFC calculation and the neg/pos merge
    have all run on top of it.
    """
    data = os.path.join(os.path.dirname(__file__), "data", "count_table.txt")
    summaries = {}
    for impl in ("cpp", "python"):
        workdir = tmp_path / impl
        workdir.mkdir()
        env = dict(os.environ, MAGECK2_RRA_IMPL=impl)
        proc = subprocess.run(
            ["mageck2", "test", "-k", data,
             "-t", "HL60.final,KBM7.final", "-c", "HL60.initial,KBM7.initial",
             "-n", "out"],
            cwd=str(workdir), env=env, capture_output=True, text=True)
        assert proc.returncode == 0, proc.stdout + proc.stderr
        summaries[impl] = (workdir / "out.gene_summary.txt").read_text()
        if impl == "python":
            assert "Python RRA implementation" in proc.stdout + proc.stderr, (
                "the switch did not take effect, so this test compared the "
                "compiled helper against itself")

    assert summaries["cpp"] == summaries["python"]
