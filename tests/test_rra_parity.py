"""Acceptance tests for the RRA C++ -> Python port.

These are the tests that decide whether the port may replace the compiled
helper. They compare an implementation against goldens recorded from the C++
(see ``tests/parity/``), gene by gene and column by column.

The file is written before the port exists, on purpose. Until
``mageck2.rra`` is importable, ``test_python_port_*`` skips; the rest of the
file still runs and does real work:

* ``test_cpp_still_reproduces_its_goldens`` re-runs the compiled helper against
  the stored goldens, so the goldens are continuously validated for as long as
  the C++ is in the tree. If the two ever diverge, the goldens were recorded
  from a different build and must be re-recorded deliberately rather than
  discovered to be stale halfway through the port.
* ``test_comparator_*`` mutation-test the comparator itself. A differential
  harness that cannot fail is worse than no harness, because it converts "we
  did not check" into "we checked and it passed". Each of these feeds the
  comparator a known-wrong candidate and requires it to object.
"""

import json
import pathlib
import shutil
import tempfile

import numpy
import pytest

from mageck2 import rra
from parity import compare, corpus, runners

GOLDENS = pathlib.Path(__file__).parent / "parity" / "goldens"
MANIFEST = GOLDENS / "manifest.json"

pytestmark = pytest.mark.skipif(
    not MANIFEST.exists(),
    reason="no goldens recorded; run `python tests/parity/record.py` while "
           "the compiled RRA helper still exists",
)

PARITY_IDS = [c.name for c in corpus.PARITY_CASES]


def manifest():
    return json.loads(MANIFEST.read_text())


def golden_text(name):
    return (GOLDENS / name / "expected.txt").read_text()


# --------------------------------------------------------------------------
# The goldens themselves
# --------------------------------------------------------------------------

@pytest.mark.parametrize("case", corpus.PARITY_CASES, ids=PARITY_IDS)
def test_every_parity_case_has_a_golden(case):
    """A case with no golden is a path the port is not actually checked on."""
    assert (GOLDENS / case.name / "expected.txt").exists(), (
        "case %r has no recorded golden" % case.name)
    assert (GOLDENS / case.name / "input.txt").exists(), (
        "case %r has no recorded input; the golden cannot be reproduced "
        "without it" % case.name)
    assert case.name in manifest()["cases"]


@pytest.mark.parametrize("case", corpus.CONTRACT_CASES,
                         ids=[c.name for c in corpus.CONTRACT_CASES])
def test_contract_cases_have_no_golden(case):
    """Contract cases must never acquire a golden.

    These reach undefined behaviour in the C++ (``control_empty`` reads
    ``control_prob_array[-1]``). Recording a golden would promote one
    build's accident to a requirement the port has to reproduce.
    """
    assert not (GOLDENS / case.name / "expected.txt").exists(), (
        "case %r is a contract case but has a golden; delete it and assert "
        "the intended behaviour instead" % case.name)


def test_manifest_records_what_produced_the_goldens():
    """Provenance, so a golden can be traced after the C++ is deleted."""
    oracle = manifest()["oracle"]
    assert oracle["fingerprint"], "no binary fingerprint recorded"
    assert "sha256:" in oracle["fingerprint"]
    assert oracle["git_commit"], "no git commit recorded"


#: Planted hits must rank in the top few percent of genes by p value. Under
#: the null their mean normalized rank would be 0.5; across the recorded
#: corpus the worst case is 0.020, so 0.05 leaves a factor of two of headroom
#: while still being ten times better than noise. Measured, not guessed --
#: a fractional "how many were called at FDR<=0.1" threshold was tried first
#: and rejected as too brittle, because with six planted hits per case it
#: quantizes to sixths and sits one gene away from failing.
MAX_MEAN_HIT_RANK = 0.05


@pytest.mark.parametrize("case", corpus.PARITY_CASES, ids=PARITY_IDS)
def test_goldens_rank_the_planted_hits_near_the_top(case):
    """Sanity on the corpus: the oracle finds the effects we planted.

    Without this, a case could be comparing two implementations that agree
    only because both return noise. The check is rank-based rather than a
    count of genes under an FDR cutoff, so it stays meaningful in cases where
    a planted hit legitimately cannot reach significance -- a gene with a
    single sgRNA has its p value bounded by one uniform percentile, so it will
    rank well without ever clearing FDR 0.1.
    """
    entry = manifest()["cases"][case.name]
    planted = set(entry["planted_hits"])
    rows, _ = compare.parse(golden_text(case.name))
    present = sorted(planted & set(rows))
    assert present, "no planted hit survived into the output"

    genes = sorted(rows)
    ranks = compare.normalized_ranks([rows[g]["p"] for g in genes])
    index = {g: i for i, g in enumerate(genes)}
    mean_rank = float(sum(ranks[index[g]] for g in present) / len(present))
    assert mean_rank <= MAX_MEAN_HIT_RANK, (
        "planted hits in case %r have mean normalized p rank %.4f; expected "
        "<= %.2f (0.5 would be indistinguishable from noise)"
        % (case.name, mean_rank, MAX_MEAN_HIT_RANK))


# --------------------------------------------------------------------------
# The oracle, while it exists
# --------------------------------------------------------------------------

@pytest.mark.skipif(shutil.which("RRA") is None,
                    reason="the compiled RRA helper is not on PATH")
@pytest.mark.parametrize("case", corpus.PARITY_CASES, ids=PARITY_IDS)
def test_cpp_still_reproduces_its_goldens(case, tmp_path):
    """The stored goldens still match what the C++ produces today.

    This keeps the goldens honest while the C++ remains in the tree, and is
    also the determinism check: RRA seeds with a fixed ``PlantSeeds(123456)``,
    so a failure here means either the helper changed or it is no longer
    reproducible -- both of which must be understood before the port lands.
    """
    corpus.build(case, tmp_path)
    produced, _ = runners.CppRunner().run(case, tmp_path)
    report = compare.compare(case.name, golden_text(case.name), produced)
    assert report.ok, report.render()


# --------------------------------------------------------------------------
# The port
# --------------------------------------------------------------------------

@pytest.mark.parametrize("case", corpus.PARITY_CASES, ids=PARITY_IDS)
def test_python_port_matches_goldens(case, tmp_path):
    """The acceptance test. Until the port exists, this skips."""
    runner = runners.PythonRunner()
    if not runner.available():
        pytest.skip("mageck2.rra does not exist yet")
    corpus.build(case, tmp_path)
    produced, _ = runner.run(case, tmp_path)
    report = compare.compare(case.name, golden_text(case.name), produced)
    assert report.ok, report.render()


@pytest.mark.parametrize("case", corpus.CONTRACT_CASES,
                         ids=[c.name for c in corpus.CONTRACT_CASES])
def test_python_port_honours_the_control_contract(case, tmp_path):
    """``--control`` matching nothing must fail loudly, not read out of bounds.

    The C++ computes ``rand_ctl_index = -1`` and reads before the start of
    ``control_prob_array`` (``RRA.cpp:761-762``, issue #51). A silent fallback
    to the uniform null would be worse than the crash, because the reported p
    values would not be the control-based ones the user asked for and nothing
    would say so. The contract is therefore: refuse to run.
    """
    runner = runners.PythonRunner()
    if not runner.available():
        pytest.skip("mageck2.rra does not exist yet")
    corpus.build(case, tmp_path)
    with pytest.raises(Exception) as excinfo:
        runner.run(case, tmp_path)
    assert "control" in str(excinfo.value).lower(), (
        "the failure did not mention the control list: %s" % excinfo.value)


# --------------------------------------------------------------------------
# The comparator, so it cannot go vacuous
# --------------------------------------------------------------------------

def _rewrite(text, row_index, column_index, value):
    """Return ``text`` with one field of one data row replaced."""
    lines = text.splitlines()
    fields = [f.strip() for f in lines[1 + row_index].split("\t")]
    fields[column_index] = value
    lines[1 + row_index] = "\t".join(fields)
    return "\n".join(lines) + "\n"


@pytest.fixture
def gold():
    """A tie-heavy golden: the case most likely to expose a bad comparison."""
    return golden_text("tie_heavy")


def test_comparator_accepts_an_identical_candidate(gold):
    report = compare.compare("tie_heavy", gold, gold)
    assert report.ok, report.render()
    assert report.deviations["p"].max_abs == 0.0
    assert report.rank_inversions == 0


def test_comparator_tolerates_differences_below_printed_precision(gold):
    """Agreement to the last printed digit must pass -- but be reported.

    ``%10.4e`` carries five significant figures, so requiring more than that
    would make the harness fail on arithmetic the file cannot even represent.
    The deviation is still surfaced so that drift inside the tolerance is
    visible rather than invisible.
    """
    rows, _ = compare.parse(gold)
    nudged = gold
    for i in range(len(rows)):
        fields = [f.strip() for f in gold.splitlines()[1 + i].split("\t")]
        nudged = _rewrite(nudged, i, 3,
                          "%10.4e" % (float(fields[3]) * (1 + 3e-5)))
    report = compare.compare("tie_heavy", gold, nudged)
    assert report.ok, report.render()
    assert report.deviations["p"].max_rel > 0.0, (
        "a real difference was reported as exactly zero")


def test_comparator_rejects_a_difference_above_printed_precision(gold):
    rows, _ = compare.parse(gold)
    first_p = float([f.strip() for f in gold.splitlines()[1].split("\t")][3])
    broken = _rewrite(gold, 0, 3, "%10.4e" % (first_p * 1.001))
    report = compare.compare("tie_heavy", gold, broken)
    assert not report.ok
    assert report.deviations["p"].n_outside == 1
    assert "OUTSIDE TOLERANCE" in report.render()


def test_comparator_rejects_a_changed_sgrna_count(gold):
    fields = [f.strip() for f in gold.splitlines()[1].split("\t")]
    broken = _rewrite(gold, 0, 5, str(int(fields[5]) + 1))
    report = compare.compare("tie_heavy", gold, broken)
    assert not report.ok
    assert report.structural, "a goodsgrna change was not reported"


def test_comparator_rejects_a_missing_gene(gold):
    lines = gold.splitlines()
    broken = "\n".join(lines[:3] + lines[4:]) + "\n"
    report = compare.compare("tie_heavy", gold, broken)
    assert not report.ok
    assert len(report.missing) == 1


def test_comparator_rejects_an_fdr_turned_into_na(gold):
    broken = _rewrite(gold, 0, 4, "NA")
    report = compare.compare("tie_heavy", gold, broken)
    assert not report.ok
    assert report.na_only_candidate, "an FDR that became NA was not reported"


def test_comparator_rejects_a_different_case_entirely(gold):
    """The crudest check: two different screens must not compare equal."""
    report = compare.compare("tie_heavy", gold, golden_text("homogeneous_5"))
    assert not report.ok
    assert report.rank_inversions > 0


def test_comparator_ignores_row_order_but_reports_it(gold):
    """Rows are sorted on a value that ties, so order is not a requirement.

    Reordering rows changes nothing scientifically -- the comparison is keyed
    by gene -- so it must still pass, while the report says the order differs
    so a real sorting change is not invisible.
    """
    lines = gold.splitlines()
    reordered = "\n".join([lines[0]] + list(reversed(lines[1:]))) + "\n"
    report = compare.compare("tie_heavy", gold, reordered)
    assert report.ok, report.render()
    assert not report.order_matches


def test_comparator_counts_a_swapped_pair_as_an_inversion(gold):
    """A p value swap between two untied genes is a real ordering change."""
    rows, order = compare.parse(gold)
    a, b = order[0], order[-1]
    assert rows[a]["p"] != rows[b]["p"], "fixture genes must not tie"
    swapped = gold
    swapped = _rewrite(swapped, 0, 3, "%10.4e" % rows[b]["p"])
    swapped = _rewrite(swapped, len(order) - 1, 3, "%10.4e" % rows[a]["p"])
    report = compare.compare("tie_heavy", gold, swapped)
    assert report.rank_inversions > 0
    assert not report.ok


# --------------------------------------------------------------------------
# The PythonRunner plumbing, before there is a port to put behind it
# --------------------------------------------------------------------------

class _StubResult:
    def __init__(self, text):
        self._text = text

    def to_text(self):
        return self._text


class _StubRraModule:
    """Minimal stand-in for the module the port will provide.

    Its only purpose is to execute ``PythonRunner.run`` so that the interface
    the port has to satisfy -- ``rank_from_file(...) -> object with
    .to_text()`` -- is exercised now rather than discovered to be wrong on the
    day the port first runs.
    """

    def __init__(self, text, recorder=None):
        self._text = text
        self.calls = [] if recorder is None else recorder

    def rank_from_file(self, path, alpha=None, control=None, flags=None):
        self.calls.append({"path": path, "alpha": alpha,
                           "control": control, "flags": flags})
        return _StubResult(self._text)


@pytest.fixture
def stub_port(monkeypatch):
    """Install a stub ``mageck2.rra`` and hand back a factory for it."""
    import sys
    import types

    def install(text):
        module = _StubRraModule(text)
        package = sys.modules.get("mageck2")
        if package is None:
            package = types.ModuleType("mageck2")
            package.__version__ = "stub"
            monkeypatch.setitem(sys.modules, "mageck2", package)
        monkeypatch.setattr(package, "rra", module, raising=False)
        monkeypatch.setitem(sys.modules, "mageck2.rra", module)
        return module

    return install


def test_python_runner_passes_the_case_through_to_the_port(stub_port, tmp_path):
    """The runner must hand the port the alpha, flags and control list.

    A runner that quietly dropped ``--control`` or the alpha would compare the
    port's default behaviour against a golden recorded with those settings,
    and the comparison would be meaningless while appearing to pass.
    """
    case = corpus.BY_NAME["control_heterogeneous"]
    module = stub_port(golden_text(case.name))
    corpus.build(case, tmp_path)

    produced, _ = runners.PythonRunner().run(case, tmp_path)

    assert len(module.calls) == 1
    call = module.calls[0]
    assert call["alpha"] == case.alpha
    assert call["flags"] == list(case.flags)
    assert call["control"], "the control list was not passed to the port"
    assert call["path"].endswith("input.txt")
    report = compare.compare(case.name, golden_text(case.name), produced)
    assert report.ok, report.render()


def test_python_runner_reports_a_wrong_port_as_a_failure(stub_port, tmp_path):
    """End-to-end proof the acceptance path can fail.

    The stub returns another case's output. If this passed, every
    ``test_python_port_matches_goldens`` result would be worthless.
    """
    case = corpus.BY_NAME["tie_heavy"]
    stub_port(golden_text("homogeneous_5"))
    corpus.build(case, tmp_path)

    produced, _ = runners.PythonRunner().run(case, tmp_path)
    report = compare.compare(case.name, golden_text(case.name), produced)
    assert not report.ok, (
        "the harness accepted a port that returned a different screen's "
        "results entirely")


# --------------------------------------------------------------------------
# The input contract: what the group column means, and what the two
# optional columns mean. None of this can be a recorded golden -- the
# long-name case makes the C++ print two rows under one label, which a
# name-keyed golden cannot align, and the optional columns are refused
# rather than reproduced.
# --------------------------------------------------------------------------

def _flat_input(rows):
    """``(items, groups, values, lists)`` from ``(sgrna, group, value)``."""
    return ([r[0] for r in rows], [r[1] for r in rows],
            [r[2] for r in rows], [0] * len(rows))


def test_group_column_is_a_list_not_a_name():
    """``P1,P2`` makes the row a member of P1 and of P2.

    The group field is comma-separated and the C++ loops over the parts
    (``fileio.cpp:107-141``, ``:266-287``). mageck2's own pathway mode relies
    on it -- ``pathwayFunc.py:155-162`` joins every pathway containing a gene
    into one field -- so reading the field as a single name invents a pathway
    and shorts each real one a member. ``pathway_memberships`` covers this
    against the helper; this states the rule in one place.
    """
    items, groups, values, lists = _flat_input([
        ("sg0", "P1,P2", 0.1), ("sg1", "P1", 0.2), ("sg2", "P2", 0.3),
        ("sg3", "P3", 0.4), ("sg4", "P3", 0.5),
    ])
    result = rra.rank(items, groups, values, lists=lists, permutation=10,
                      permutation_specified=True)
    counts = dict(zip(result.names, result.items.tolist()))
    assert counts == {"P1": 2, "P2": 2, "P3": 2}, (
        "expected sg0 to count in both P1 and P2; got %r" % (counts,))


def test_a_repeated_group_name_is_not_deduplicated():
    """``P1,P1`` gives P1 the row twice.

    The C++ increments ``itemNum`` once per occurrence and appends the item
    once per occurrence, so the count exceeds the number of distinct rows.
    Verified against the helper: on a probe input, ``PW1,PW2`` plus
    ``PW1,PW1`` plus ``PW2`` gave PW1 three items and PW2 two.
    """
    items, groups, values, lists = _flat_input([
        ("sg0", "P1,P1", 0.1), ("sg1", "P1", 0.2),
    ])
    result = rra.rank(items, groups, values, lists=lists, permutation=10,
                      permutation_specified=True)
    counts = dict(zip(result.names, result.items.tolist()))
    assert counts == {"P1": 3}, (
        "a repeated sub-name must not be collapsed; got %r" % (counts,))


def test_grouping_keys_on_the_full_name_not_the_stored_label():
    """Two names agreeing to ``MAX_NAME_LEN - 1`` stay two groups.

    ``copyName`` (``fileio.cpp:44-58``) truncates the label it stores and says
    so in its warning -- "Grouping is unaffected -- it keys on the full name"
    -- while ``groupNames`` is keyed on the untruncated string. Keying on the
    truncated form instead merges distinct genes, which changes their scores
    and the size of the population the permutation and FDR are computed over.

    Measured against the helper on a probe input built exactly this way: it
    printed two rows under the same 511-character label, with 4 and 6 sgRNAs
    and lo-values 3.2594e-03 and 9.1451e-01. The port reproduced both
    numbers; keying on the truncated name produced one row of 10 instead.

    The row count alone does not catch it, which is worth knowing: in that
    probe the merge lost one row and the unsplit comma fields gained one, so
    both bugs together left the total unchanged at 24.
    """
    prefix = "x" * (rra.MAX_NAME_LEN - 1)
    long_a, long_b = prefix + "_AAA", prefix + "_BBB"
    assert long_a[:rra.MAX_NAME_LEN - 1] == long_b[:rra.MAX_NAME_LEN - 1]

    rows = [("sgA%d" % i, long_a, 0.1 * (i + 1)) for i in range(4)]
    rows += [("sgB%d" % i, long_b, 0.5 + 0.05 * i) for i in range(6)]
    items, groups, values, lists = _flat_input(rows)
    result = rra.rank(items, groups, values, lists=lists, permutation=10,
                      permutation_specified=True)

    assert len(result.names) == 2, (
        "the two names must stay two groups; got %r group(s)"
        % (len(result.names),))
    assert sorted(result.items.tolist()) == [4, 6]
    # Both labels are the truncated form, which is what the C++ writes: the
    # collision is in the output, not in the grouping.
    assert set(result.names) == {prefix}
    assert all(len(n) == rra.MAX_NAME_LEN - 1 for n in result.names)


def _write_input(path, rows, header="sgrna\tsymbol\tpool\tvalue"):
    path.write_text(header + "\n" + "".join(r + "\n" for r in rows))
    return path


def test_read_input_accepts_the_neutral_optional_columns():
    """The six-column form mageck2 itself writes is accepted unchanged.

    ``crisprFunction.py:505`` and ``:529`` write ``prob`` as the literal 1 and
    only emit rows whose ``chosen`` is 1, so the shipped pipeline must keep
    working. If this test ever fails, the rejection below has become a
    regression in the mainline path rather than a guard on an unsupported one.
    """
    with tempfile.TemporaryDirectory() as tmp:
        path = _write_input(
            pathlib.Path(tmp) / "input.txt",
            ["sg0\tG1\tlist\t0.1\t1\t1", "sg1\tG1\tlist\t0.2\t1\t1"],
            header="sgrna\tsymbol\tpool\tp.low\tprob\tchosen")
        items, groups, values, pools = rra.read_input(str(path))
    assert items == ["sg0", "sg1"]
    assert groups == ["G1", "G1"]
    assert list(values) == [0.1, 0.2]


@pytest.mark.parametrize("row,needle", [
    ("sg0\tG1\tlist\t0.1\t0.5\t1", "probabilities other than 1"),
    ("sg0\tG1\tlist\t0.1\t1\t0", "unchosen sgRNAs"),
])
def test_read_input_refuses_what_it_cannot_reproduce(row, needle):
    """A meaningful ``prob`` or ``chosen`` is refused, not ignored.

    Both columns change the answer in the C++: ``prob`` below 1 switches to
    ``ComputeLoValue_Prob`` (``RRA.cpp:525``), and ``chosen`` of 0 drops the
    row from its pool's percentiles and from its gene's lo-value while still
    counting it in that gene's sgRNA total. Reading past them returns
    different statistics for an accepted input with nothing said about it,
    which is the one outcome a port must not have. Refusing is a deliberate
    narrowing of what the Python backend claims to support.
    """
    with tempfile.TemporaryDirectory() as tmp:
        path = _write_input(pathlib.Path(tmp) / "input.txt", [row])
        with pytest.raises(ValueError) as excinfo:
            rra.read_input(str(path))
    assert needle in str(excinfo.value)
    assert "MAGECK2_RRA_IMPL=cpp" in str(excinfo.value), (
        "the error must say how to get the supported behaviour")


def test_null_chunking_is_bounded_by_elements_not_rows():
    """A wide gene must not defeat the chunk limit.

    The bound exists to cap peak memory while the null is built. Counting rows
    instead of elements makes the cap depend on gene width: at 2**18 rows a
    100-sgRNA gene would materialize a 262,144 x 100 float64 array, 210 MB,
    and the reduction holds several arrays of that shape at once. Tiling and
    saturation libraries really do put that many sgRNAs on one gene.
    """
    chunk = 1 << 12
    widths = []

    class _Recorder:
        """Counts the draws asked for, which is the chunk's element count."""

        def uniform(self, n):
            widths.append(n)
            return numpy.full(n, 0.5, dtype=float)

    rra._null_lo_values(_Recorder(), rows=500, k=100, alpha=0.25,
                        control_values=None, chunk=chunk)
    assert widths, "no draws were requested"
    assert max(widths) <= chunk, (
        "a chunk asked for %d elements with a limit of %d"
        % (max(widths), chunk))
    # And it still has to cover every row exactly once.
    assert sum(widths) == 500 * 100
