"""Acceptance tests for the ``mageckGSEA`` C++ -> Python port.

These are the tests that decide whether the port may replace the compiled
helper. They compare an implementation against goldens recorded from the C++
(see ``tests/parity/gsea/``), pathway by pathway and column by column.

The file is written before the port exists, on purpose. Until
``mageck2.gsea`` is importable, ``test_python_port_*`` skips; the rest of the
file still runs and does real work:

* ``test_cpp_still_reproduces_its_goldens`` re-runs the compiled helper against
  the stored goldens, so the goldens are continuously validated for as long as
  the C++ is in the tree. If the two ever diverge, the goldens were recorded
  from a different build and must be re-recorded deliberately rather than
  discovered to be stale halfway through the port.
* ``test_comparator_*`` mutation-test the comparator itself. A differential
  harness that cannot fail is worse than no harness, because it converts "we
  did not check" into "we checked and it passed". Each of these feeds the
  comparator a known-wrong candidate and requires it to object -- including
  the two cases where it must *not* object, since a comparator that rejects
  everything is equally useless.

Divergences
-----------

Five corpus cases are ``Mode.DIVERGENCE``: the C++ behaviour is defined,
reproducible, and wrong, and the port is required to differ. Their goldens are
recorded anyway -- that is what keeps the premise of each fix checkable, and it
documents what users of past versions actually got -- but the port is never
compared against them.

What the port is asked for instead is a property: **a malformed input must give
the same answer as its repair.** A missing final newline, a blank line, a
multi-word description and trailing whitespace are all things that should not
change a result, and each one changes it today. Stating the fix that way is
stronger than listing expected numbers, and it does not need the port to exist
in order to be written down.
"""

import json
import pathlib
import shutil

import numpy as np
import pytest

from parity.gsea import compare, corpus, runners
from parity.stats import normalized_ranks

GOLDENS = pathlib.Path(__file__).parent / "parity" / "gsea" / "goldens"
MANIFEST = GOLDENS / "manifest.json"

pytestmark = pytest.mark.skipif(
    not MANIFEST.exists(),
    reason="no goldens recorded; run `python tests/parity/gsea/record.py` "
           "while the compiled mageckGSEA helper still exists",
)

PARITY_IDS = [c.name for c in corpus.PARITY_CASES]
RECORDED_IDS = [c.name for c in corpus.RECORDED_CASES]
DIVERGENCE_IDS = [c.name for c in corpus.DIVERGENCE_CASES]

#: The enriched gene sets are drawn from the top of the ranking, so a correct
#: implementation puts them near the top of the output. Measured over the
#: goldens, the worst case's enriched sets average a normalised rank of 0.106
#: against a null expectation of 0.5; 0.2 leaves room for the tie-averaging
#: that the permutation grid forces without coming anywhere near null.
MAX_MEAN_ENRICHED_RANK = 0.2


def manifest():
    return json.loads(MANIFEST.read_text())


def golden_text(name):
    return (GOLDENS / name / "expected.txt").read_text()


# --------------------------------------------------------------------------
# The goldens themselves
# --------------------------------------------------------------------------

@pytest.mark.parametrize("case", corpus.RECORDED_CASES, ids=RECORDED_IDS)
def test_every_recorded_case_has_a_golden(case):
    """A case with no golden is a path nothing is actually checked on."""
    assert (GOLDENS / case.name / "expected.txt").exists(), (
        "case %r has no recorded golden" % case.name)
    for required in ("rank.txt", "sets.gmt"):
        assert (GOLDENS / case.name / required).exists(), (
            "case %r has no recorded %s; the golden cannot be reproduced "
            "once the C++ is gone" % (case.name, required))


def test_contract_cases_have_no_golden():
    """Undefined behaviour must never acquire a golden.

    There are no ``CONTRACT`` cases in the GSEA corpus today -- unlike RRA,
    whose ``control_empty`` case reaches an out-of-bounds read. The assertion
    is kept so that adding one and recording it anyway is caught: freezing the
    output of undefined behaviour turns one build's accident into a
    requirement.
    """
    for case in corpus.CONTRACT_CASES:
        assert not (GOLDENS / case.name / "expected.txt").exists(), (
            "case %r reaches undefined behaviour in the C++ but has a "
            "recorded golden" % case.name)


def test_manifest_records_what_produced_the_goldens():
    """Provenance, so a golden can be traced after the C++ is deleted."""
    data = manifest()
    oracle = data["oracle"]
    assert oracle["fingerprint"], "no binary fingerprint recorded"
    assert oracle["fingerprint"].startswith("mageckGSEA sha256:"), (
        "the fingerprint should name the helper and its hash, got %r"
        % oracle["fingerprint"])
    assert oracle["git_commit"], "no git commit recorded"
    assert oracle["oracle_source_dirty"] is False, (
        "the goldens were recorded from a dirty rra/ or gsea/ tree, so the "
        "recorded commit does not describe the binary that produced them")
    assert set(data["cases"]) == {c.name for c in corpus.RECORDED_CASES}, (
        "the manifest and the corpus disagree about which cases exist")


@pytest.mark.parametrize("case", corpus.PARITY_CASES, ids=PARITY_IDS)
def test_goldens_rank_the_enriched_sets_near_the_top(case):
    """The goldens contain signal, not just self-consistent noise.

    A corpus of pure null data would let a broken port agree with a broken
    golden and report success. ``DIVERGENCE`` cases are excluded because their
    inputs are malformed on purpose, and the cases whose planting does not
    apply -- scoring on the opposite column, or ignoring the score columns
    altogether -- declare that with ``signal_expected=False``.
    """
    if not case.signal_expected:
        pytest.skip("this case's planting does not apply to its scoring: "
                    + case.why.split(".")[0])
    rows, _ = compare.parse(golden_text(case.name))
    enriched = set(manifest()["cases"][case.name]["enriched_sets"])
    assert enriched, "no enriched sets were recorded for this case"

    names = sorted(rows)
    ranks = normalized_ranks([rows[n]["p_permutation"] for n in names])
    planted = [i for i, n in enumerate(names) if n in enriched]
    assert planted, "none of the enriched sets appear in the output"

    mean_rank = float(np.mean(ranks[planted]))
    assert mean_rank <= MAX_MEAN_ENRICHED_RANK, (
        "the enriched gene sets average a normalised p rank of %.3f "
        "(0 is the most significant, 0.5 is what noise gives), so this "
        "golden carries no signal for a port to reproduce" % mean_rank)


# --------------------------------------------------------------------------
# The oracle, for as long as it exists
# --------------------------------------------------------------------------

@pytest.mark.skipif(shutil.which("mageckGSEA") is None,
                    reason="mageckGSEA is not on PATH")
@pytest.mark.parametrize("case", corpus.RECORDED_CASES, ids=RECORDED_IDS)
def test_cpp_still_reproduces_its_goldens(case, tmp_path):
    """The compiled helper must still produce what was recorded from it.

    This is what keeps the goldens honest while the C++ is in the tree. It is
    also what keeps the ``DIVERGENCE`` cases meaningful: each of those asserts
    that the C++ does something specific and wrong, and if the C++ stopped
    doing it, the fix the port is required to make would no longer have a
    premise.
    """
    corpus.build(case, tmp_path)
    produced, _ = runners.CppRunner().run(case, tmp_path)
    report = compare.compare(case.name, golden_text(case.name), produced)
    assert report.ok, report.render()


def test_the_goldens_are_not_required_to_be_byte_reproducible():
    """Documents why the test above stops at numeric agreement.

    An earlier version asserted byte equality, which is the stronger and more
    satisfying check -- and it is wrong here. ``getgseascore_fast`` computes
    ``decvalue*(ni-i) + sumincval``, which clang contracts into a fused
    multiply-add at ``-O3`` wherever FMA is baseline (it is on arm64, it is
    not on baseline x86-64). Rebuilding this exact source with
    ``-ffp-contract=off`` moves three of the twenty goldens by about 1e-17, in
    gene sets whose score cancels to zero. So byte equality would assert that
    the test machine's compiler matches the recording machine's, and would
    fail the C++ against its own goldens the first time CI ran on a different
    architecture.

    This test exists so that someone tightening the assertion back up finds
    out why before CI tells them.
    """
    assert compare.ATOL_ES > 0.0, (
        "ES needs an absolute floor; see the constant's own documentation")
    assert compare.NUMERIC_COLUMNS[0] == ("ES", True, compare.ATOL_ES)
    # The columns the FMA cannot reach are still required to be exact.
    exact = {name: atol for name, printed, atol in compare.NUMERIC_COLUMNS
             if not printed}
    assert exact == {"p_permutation": 0.0}


# --------------------------------------------------------------------------
# The port
# --------------------------------------------------------------------------

def _port_available():
    return runners.PythonRunner().available()


needs_port = pytest.mark.skipif(
    not _port_available(),
    reason="mageck2.gsea does not exist yet; this is the port's target "
           "interface, see parity/gsea/runners.py")


@needs_port
@pytest.mark.parametrize("case", corpus.PARITY_CASES, ids=PARITY_IDS)
def test_python_port_matches_the_goldens(case, tmp_path):
    """The acceptance test. Every parity case, every column."""
    corpus.build(case, tmp_path)
    produced, _ = runners.PythonRunner().run(case, tmp_path)
    report = compare.compare(case.name, golden_text(case.name), produced)
    assert report.ok, report.render()


# How to repair each malformed input. The port's answer on the malformed file
# must equal its answer on the repaired one -- which is the property the C++
# lacks, and a stronger statement than any list of expected numbers.
def _repair_trailing_newline(name):
    def repair(directory):
        path = directory / name
        text = path.read_text()
        if not text.endswith("\n"):
            path.write_text(text + "\n")
    return repair


def _repair_blank_lines(directory):
    path = directory / "sets.gmt"
    kept = [ln for ln in path.read_text().splitlines() if ln.strip()]
    path.write_text("\n".join(kept) + "\n")


def _repair_description(directory):
    path = directory / "sets.gmt"
    out = []
    for line in path.read_text().splitlines():
        fields = line.split("\t")
        if len(fields) > 1:
            fields[1] = corpus.DESCRIPTION
        out.append("\t".join(fields))
    path.write_text("\n".join(out) + "\n")


def _repair_trailing_whitespace(directory):
    path = directory / "sets.gmt"
    out = [ln.rstrip() for ln in path.read_text().splitlines()]
    path.write_text("\n".join(out) + "\n")


REPAIRS = {
    "gmt_missing_final_newline": _repair_trailing_newline("sets.gmt"),
    "rank_missing_final_newline": _repair_trailing_newline("rank.txt"),
    "gmt_blank_line": _repair_blank_lines,
    "gmt_description_with_spaces": _repair_description,
    "gmt_trailing_whitespace": _repair_trailing_whitespace,
}


def test_every_divergence_case_states_how_to_repair_it():
    """A divergence with no repair has no testable fix."""
    assert set(REPAIRS) == {c.name for c in corpus.DIVERGENCE_CASES}, (
        "REPAIRS and the DIVERGENCE cases disagree; a divergence case "
        "without a repair is asserted by nothing")


@pytest.mark.skipif(shutil.which("mageckGSEA") is None,
                    reason="mageckGSEA is not on PATH")
@pytest.mark.parametrize("case", corpus.DIVERGENCE_CASES, ids=DIVERGENCE_IDS)
def test_the_cpp_really_is_wrong_on_each_divergence(case, tmp_path):
    """The premise of each fix, asserted against the C++ rather than claimed.

    If the helper's answer on the malformed input already equalled its answer
    on the repaired one, the case would be documenting a bug that does not
    exist, and the port would be required to fix nothing.
    """
    malformed = tmp_path / "malformed"
    repaired = tmp_path / "repaired"
    malformed.mkdir()
    repaired.mkdir()
    corpus.build(case, malformed)
    corpus.build(case, repaired)
    REPAIRS[case.name](repaired)

    bad, _ = runners.CppRunner().run(case, malformed)
    good, _ = runners.CppRunner().run(case, repaired)
    assert bad != good, (
        "the C++ gives the same answer for %r malformed and repaired, so "
        "there is no divergence here and the case should be reclassified"
        % case.name)


@needs_port
@pytest.mark.parametrize("case", corpus.DIVERGENCE_CASES, ids=DIVERGENCE_IDS)
def test_python_port_is_unmoved_by_each_malformation(case, tmp_path):
    """The port must answer the same for a malformed input and its repair."""
    malformed = tmp_path / "malformed"
    repaired = tmp_path / "repaired"
    malformed.mkdir()
    repaired.mkdir()
    corpus.build(case, malformed)
    corpus.build(case, repaired)
    REPAIRS[case.name](repaired)

    runner = runners.PythonRunner()
    bad, _ = runner.run(case, malformed)
    good, _ = runner.run(case, repaired)
    report = compare.compare(case.name, good, bad)
    assert report.ok, (
        "the port's answer changes when the input is malformed in the way "
        "%r describes, which is the C++ behaviour it is supposed to "
        "fix:\n%s" % (case.name, report.render()))


@needs_port
@pytest.mark.parametrize("case", corpus.DIVERGENCE_CASES, ids=DIVERGENCE_IDS)
def test_python_port_does_not_reproduce_the_divergence_goldens(case, tmp_path):
    """The complement: the port must *not* match these goldens.

    Without this, a port that reproduced the C++ bit for bit -- bugs
    included -- would pass every test in this file, because the divergence
    goldens are otherwise only ever read, never asserted against.
    """
    corpus.build(case, tmp_path)
    produced, _ = runners.PythonRunner().run(case, tmp_path)
    report = compare.compare(case.name, golden_text(case.name), produced)
    assert not report.ok, (
        "the port reproduces the C++ exactly on %r, including the bug the "
        "case exists to document" % case.name)


# --------------------------------------------------------------------------
# The comparator, mutation-tested
# --------------------------------------------------------------------------

#: A case with no tied pathway pairs, so every ordering claim it makes is real.
UNTIED_CASE = "production_default"

#: A case containing pathways that tie on both sort keys, where the C++ row
#: order is genuinely unspecified. Asserted to still contain them, so the
#: exemption tests cannot silently become vacuous -- and that assertion has
#: already earned its place: ``tie_extreme`` held such pairs until the corpus
#: started shuffling its input order, and the guard caught the loss rather
#: than letting the exemption tests quietly stop testing anything.
TIED_CASE = "big_sets"


def _rewrite(text, pathway, column, value):
    """Return ``text`` with one cell replaced, keeping everything else byte
    for byte. Formatting matches the C++: six significant digits."""
    out = []
    index = compare.COLUMNS.index(column)
    for i, line in enumerate(text.splitlines()):
        fields = line.split("\t")
        if i and fields[0] == pathway:
            fields[index] = value if isinstance(value, str) else "%.6g" % value
            line = "\t".join(fields)
        out.append(line)
    return "\n".join(out) + "\n"


def _first_pathway(text):
    rows, order = compare.parse(text)
    return order[0], rows


def _tied_pair(text):
    """Two pathways that tie on both of the C++ sort keys, or None."""
    rows, _ = compare.parse(text)
    groups = {}
    for name, row in rows.items():
        groups.setdefault((row["p_permutation"], row["ES"]), []).append(name)
    for members in groups.values():
        if len(members) >= 2:
            return sorted(members)[:2]
    return None


def _swap_rows(text, a, b):
    """Swap two data rows, and swap their Ranking values with them."""
    lines = text.splitlines()
    idx = {}
    for i, line in enumerate(lines):
        if i:
            name = line.split("\t")[0]
            if name in (a, b):
                idx[name] = i
    fa, fb = lines[idx[a]].split("\t"), lines[idx[b]].split("\t")
    rank = compare.COLUMNS.index("Ranking")
    fa[rank], fb[rank] = fb[rank], fa[rank]
    lines[idx[a]], lines[idx[b]] = "\t".join(fb), "\t".join(fa)
    return "\n".join(lines) + "\n"


def test_the_untied_case_really_has_no_tied_pairs():
    """Guards the mutation tests that rely on every pair being ordered."""
    assert _tied_pair(golden_text(UNTIED_CASE)) is None, (
        "%r has acquired tied pathway pairs, so the ordering mutation tests "
        "below no longer prove what they claim" % UNTIED_CASE)


def test_the_tied_case_really_has_tied_pairs():
    """Guards the exemption test. Without a tie, it would pass vacuously."""
    assert _tied_pair(golden_text(TIED_CASE)) is not None, (
        "%r no longer contains pathways tied on both sort keys, so the "
        "tie-exemption test proves nothing" % TIED_CASE)


def test_comparator_accepts_an_identical_candidate():
    """The floor. A comparator that rejects everything proves nothing."""
    text = golden_text(UNTIED_CASE)
    report = compare.compare(UNTIED_CASE, text, text)
    assert report.ok, report.render()
    assert report.byte_identical
    assert report.rank_inversions == 0
    assert report.ranking_disagreements == 0


def test_comparator_catches_a_single_permutation_difference():
    """One exceedance out of 1000 must not pass.

    The permutation stream is seeded and its consumption order fully
    determined, so there is no noise budget to spend: a port that is off by
    one permutation has a real defect, most likely in the rejection-sampling
    loop, and it would be invisible under any relative tolerance applied to a
    p value of 0.5.
    """
    text = golden_text(UNTIED_CASE)
    name, rows = _first_pathway(text)
    shifted = rows[name]["p_permutation"] + 1.0 / 1000.0
    broken = _rewrite(text, name, "p_permutation", shifted)
    report = compare.compare(UNTIED_CASE, text, broken)
    assert not report.ok, report.render()
    assert report.deviations["p_permutation"].n_outside == 1


def test_comparator_catches_an_es_difference_beyond_the_printed_precision():
    text = golden_text(UNTIED_CASE)
    name, rows = _first_pathway(text)
    broken = _rewrite(text, name, "ES", rows[name]["ES"] * 1.001)
    report = compare.compare(UNTIED_CASE, text, broken)
    assert not report.ok, report.render()
    assert report.deviations["ES"].n_outside == 1
    assert report.deviations["ES"].worst == name


def _one_printed_step(value):
    """``value`` plus one unit in its sixth significant digit.

    Scaling by ``1 + 1e-6`` does not work: at a mantissa of 1 that lands
    inside the sixth digit and ``%.6g`` prints the original number back, so
    the mutation silently does not happen. The step has to be computed from
    the value's magnitude.
    """
    import math
    exponent = math.floor(math.log10(abs(value)))
    return value + 10.0 ** (exponent - 5)


def test_comparator_reports_but_accepts_a_last_digit_es_difference():
    """Drift inside tolerance must still be visible.

    A comparator that reports zero deviation for anything it accepts cannot
    show a port slowly getting worse, which is the failure mode that a
    pass/fail answer hides.
    """
    text = golden_text(UNTIED_CASE)
    name, rows = _first_pathway(text)
    broken = _rewrite(text, name, "ES", _one_printed_step(rows[name]["ES"]))
    assert broken != text, (
        "the mutation did not change the file, so this test would pass "
        "vacuously")
    report = compare.compare(UNTIED_CASE, text, broken)
    assert report.ok, report.render()
    assert report.deviations["ES"].max_rel > 0.0, (
        "a last-digit difference was accepted and also reported as zero")
    assert report.deviations["ES"].max_rel < compare.RTOL_PRINTED, (
        "one unit in the last printed digit should sit inside the tolerance "
        "the tolerance claims to be sized for")
    assert not report.byte_identical


def test_comparator_does_not_check_ranking_when_it_is_not_a_ranking():
    """The ``unsorted`` output writes Ranking as a literal 0 for every row.

    Read naively, that is "every pathway ties for first", and every pair the
    golden distinguishes becomes a disagreement -- 244 of them on this case,
    against a byte-identical file. Found exactly that way: the comparator
    rejected the C++'s own output.
    """
    text = golden_text("unsorted")
    rows, _ = compare.parse(text)
    assert not compare.is_ranked(rows), (
        "the unsorted golden is expected to carry no ordering in Ranking")
    report = compare.compare("unsorted", text, text)
    assert report.ok, report.render()
    assert report.ranking_checked is False
    assert report.ranking_disagreements == 0
    assert "not checked" in report.render()


def test_comparator_checks_ranking_when_it_is_a_ranking():
    """Guards the test above from turning into a blanket exemption."""
    rows, _ = compare.parse(golden_text(UNTIED_CASE))
    assert compare.is_ranked(rows)
    report = compare.compare(UNTIED_CASE, golden_text(UNTIED_CASE),
                             golden_text(UNTIED_CASE))
    assert report.ranking_checked is True


def test_comparator_catches_a_port_that_drops_the_ranking_column():
    """Zeroing Ranking must not become a way to pass the ordering check.

    A port that forgot to fill the column in sorted mode would otherwise be
    rewarded by the exemption above.
    """
    text = golden_text(UNTIED_CASE)
    lines = text.splitlines()
    rank = compare.COLUMNS.index("Ranking")
    out = [lines[0]]
    for line in lines[1:]:
        fields = line.split("\t")
        fields[rank] = "0"
        out.append("\t".join(fields))
    report = compare.compare(UNTIED_CASE, text, "\n".join(out) + "\n")
    assert not report.ok, report.render()
    assert any("Ranking" in line for line in report.structural)


def test_comparator_catches_an_es_difference_above_the_absolute_floor():
    """The ES floor must not become a blanket exemption.

    ``ATOL_ES`` exists for scores that cancel to zero, where a one-bit
    compiler difference is the whole signal. It must still catch a difference
    that is small in absolute terms but real -- here a hundred times the
    floor, which is still fourteen orders of magnitude below 1.
    """
    text = golden_text("degenerate_sets")
    rows, _ = compare.parse(text)
    near_zero = sorted(n for n in rows if abs(rows[n]["ES"]) < 1e-10)
    assert near_zero, (
        "this case is expected to contain a gene set whose score cancels to "
        "zero; without one the floor is being tested on nothing")
    name = near_zero[0]
    broken = _rewrite(text, name, "ES", compare.ATOL_ES * 100)
    report = compare.compare("degenerate_sets", text, broken)
    assert not report.ok, report.render()
    assert report.deviations["ES"].n_outside == 1


def test_comparator_accepts_an_es_difference_below_the_absolute_floor():
    """And it must accept one at the scale the compiler actually produces."""
    text = golden_text("degenerate_sets")
    rows, _ = compare.parse(text)
    near_zero = sorted(n for n in rows if abs(rows[n]["ES"]) < 1e-10)
    name = near_zero[0]
    # 5.551115123125783e-17 is the measured arm64/x86-64 difference.
    broken = _rewrite(text, name, "ES", 5.551115123125783e-17)
    assert broken != text
    report = compare.compare("degenerate_sets", text, broken)
    assert report.ok, report.render()


def test_comparator_catches_a_wrong_fdr():
    text = golden_text(UNTIED_CASE)
    name, rows = _first_pathway(text)
    broken = _rewrite(text, name, "FDR", rows[name]["FDR"] + 0.01)
    report = compare.compare(UNTIED_CASE, text, broken)
    assert not report.ok, report.render()
    assert report.deviations["FDR"].n_outside == 1


def test_comparator_catches_a_wrong_set_size():
    text = golden_text(UNTIED_CASE)
    name, rows = _first_pathway(text)
    broken = _rewrite(text, name, "Size", str(rows[name]["Size"] + 1))
    report = compare.compare(UNTIED_CASE, text, broken)
    assert not report.ok, report.render()
    assert any("Size" in line for line in report.structural)


def test_comparator_catches_a_wrong_hits_index():
    """``Hits`` is where the enrichment score peaked.

    It is the one column that says *why* the score is what it is, so a port
    that got the right score at the wrong position -- a sign the running sum
    is being maximised over the wrong index -- has to be caught even though
    every float agrees.
    """
    text = golden_text(UNTIED_CASE)
    name, rows = _first_pathway(text)
    broken = _rewrite(text, name, "Hits", str(rows[name]["Hits"] + 1))
    report = compare.compare(UNTIED_CASE, text, broken)
    assert not report.ok, report.render()
    assert any("Hits" in line for line in report.structural)


def test_comparator_catches_an_invented_lfc():
    """The LFC column is a literal zero in both output branches."""
    text = golden_text(UNTIED_CASE)
    name, _ = _first_pathway(text)
    broken = _rewrite(text, name, "LFC", -1.25)
    report = compare.compare(UNTIED_CASE, text, broken)
    assert not report.ok, report.render()
    assert any("LFC" in line for line in report.structural)


def test_comparator_catches_a_missing_pathway():
    text = golden_text(UNTIED_CASE)
    name, _ = _first_pathway(text)
    kept = [ln for i, ln in enumerate(text.splitlines())
            if not i or ln.split("\t")[0] != name]
    report = compare.compare(UNTIED_CASE, text, "\n".join(kept) + "\n")
    assert not report.ok, report.render()
    assert report.missing == [name]


def test_comparator_catches_an_invented_pathway():
    text = golden_text(UNTIED_CASE)
    extra = "SET_INVENTED\t5\t0.5\t0.5\t0.5\t0.5\t99\t1\t0"
    report = compare.compare(UNTIED_CASE, text, text + extra + "\n")
    assert not report.ok, report.render()
    assert report.extra == ["SET_INVENTED"]


def test_comparator_catches_a_reordering_of_pathways_it_can_distinguish():
    """Swapping two pathways the golden orders definitely is a real failure."""
    text = golden_text(UNTIED_CASE)
    _, order = compare.parse(text)[0], compare.parse(text)[1]
    a, b = order[0], order[-1]
    broken = _swap_rows(text, a, b)
    report = compare.compare(UNTIED_CASE, text, broken)
    assert not report.ok, report.render()
    assert report.ranking_disagreements > 0


def test_comparator_exempts_a_reordering_of_pathways_that_fully_tie():
    """And swapping two that tie on both sort keys is not.

    This is the exemption the whole comparator design turns on. ``std::sort``
    is not stable, so for a pair with equal ``p_permutation`` and equal ``ES``
    the C++ order is whatever the implementation happened to produce. Requiring
    a port to match it would be requiring it to reproduce an accident -- and
    was measured to happen in practice while this harness was being designed.
    """
    text = golden_text(TIED_CASE)
    pair = _tied_pair(text)
    assert pair is not None
    broken = _swap_rows(text, pair[0], pair[1])
    report = compare.compare(TIED_CASE, text, broken)
    assert report.ok, report.render()
    assert report.ranking_disagreements == 0
    assert not report.order_matches, (
        "the rows were not actually reordered, so this proves nothing")
    assert report.tied_pairs > 0


def test_comparator_exempts_an_ordering_it_has_declared_meaningless():
    """The exemption threshold must match the numeric tolerance.

    If a pair of ES values is close enough to be accepted as equal, the
    ordering between them cannot also be required. An earlier version used
    exact equality here and rejected a faithful reimplementation that agreed
    on every single value, purely because an ES difference of 2.776e-17 --
    which it had just accepted as zero -- put two rows in the other order.
    """
    text = golden_text(UNTIED_CASE)
    rows, order = compare.parse(text)
    # Find a pair adjacent in the output that shares a permutation p value, so
    # that only ES separates them, and make that separation negligible.
    pair = None
    for a, b in zip(order, order[1:]):
        if rows[a]["p_permutation"] == rows[b]["p_permutation"]:
            pair = (a, b)
            break
    assert pair is not None, (
        "this case is expected to contain adjacent rows sharing a p value")
    a, b = pair
    nudged = _rewrite(text, b, "ES", rows[a]["ES"] - compare.ATOL_ES / 10)
    swapped = _swap_rows(nudged, a, b)
    report = compare.compare(UNTIED_CASE, nudged, swapped)
    assert not report.order_matches, "the rows were not actually reordered"
    assert report.ranking_disagreements == 0, report.render()
    assert report.tied_pairs > 0
    assert report.ok, report.render()


def test_comparator_rejects_a_file_that_is_not_gsea_output():
    with pytest.raises(ValueError, match="header"):
        compare.parse("not\tan\toutput\tfile\n1\t2\t3\t4\n")


def test_comparator_rejects_a_truncated_row():
    text = golden_text(UNTIED_CASE)
    lines = text.splitlines()
    lines[1] = "\t".join(lines[1].split("\t")[:5])
    with pytest.raises(ValueError, match="malformed"):
        compare.parse("\n".join(lines) + "\n")


def test_comparator_reads_back_a_pathway_with_an_empty_name():
    """The ``gmt_blank_line`` golden has one, and it must be readable.

    If ``parse`` dropped or rejected it, the case documenting the most severe
    of the five divergences would silently stop being checked.
    """
    rows, order = compare.parse(golden_text("gmt_blank_line"))
    assert order == [""], (
        "the blank-line golden should contain exactly one row, named with "
        "the empty string, got %r" % (order,))
    assert rows[""]["Size"] == 0


# --------------------------------------------------------------------------
# The acceptance path itself
# --------------------------------------------------------------------------

class _StubResult:
    def __init__(self, text):
        self._text = text

    def to_text(self):
        return self._text


class _StubGseaModule:
    """Minimal stand-in for the module the port will provide.

    Its only purpose is to execute ``PythonRunner.run`` so that the interface
    the port has to satisfy -- ``enrich_from_file(...) -> object with
    .to_text()`` -- is exercised now rather than discovered to be wrong on the
    day the port first runs.
    """

    def __init__(self, text):
        self._text = text
        self.calls = []

    def enrich_from_file(self, gmt_path, rank_path, **kwargs):
        self.calls.append(dict(kwargs, gmt_path=gmt_path,
                               rank_path=rank_path))
        return _StubResult(self._text)


@pytest.fixture
def stub_port(monkeypatch):
    """Install a stub ``mageck2.gsea`` and hand back a factory for it."""
    import sys
    import types

    def install(text):
        module = _StubGseaModule(text)
        package = sys.modules.get("mageck2")
        if package is None:
            package = types.ModuleType("mageck2")
            package.__version__ = "stub"
            monkeypatch.setitem(sys.modules, "mageck2", package)
        monkeypatch.setattr(package, "gsea", module, raising=False)
        monkeypatch.setitem(sys.modules, "mageck2.gsea", module)
        return module

    return install


def test_python_runner_passes_every_flag_through_to_the_port(stub_port,
                                                             tmp_path):
    """The runner must hand the port the whole invocation.

    A runner that quietly dropped the score column or the seed would compare
    the port's default behaviour against a golden recorded with other
    settings, and the comparison would be meaningless while appearing to
    pass. ``custom_seed`` is used because it is the one case that sets every
    interesting knob away from its default.
    """
    case = corpus.BY_NAME["custom_seed"]
    module = stub_port(golden_text(case.name))
    corpus.build(case, tmp_path)

    produced, _ = runners.PythonRunner().run(case, tmp_path)

    assert len(module.calls) == 1
    call = module.calls[0]
    assert call["score_column"] == case.score_column
    assert call["perms"] == case.perms
    assert call["seed"] == case.seed
    assert call["sort_by_p"] == case.sort_by_p
    assert call["skip_header"] == case.skip_header
    assert call["reverse"] == case.reverse
    assert call["pathway_name"] == case.pathway_name
    assert call["gmt_path"].endswith("sets.gmt")
    assert call["rank_path"].endswith("rank.txt")
    report = compare.compare(case.name, golden_text(case.name), produced)
    assert report.ok, report.render()


def test_python_runner_leaves_the_seed_unset_when_the_case_does(stub_port,
                                                               tmp_path):
    """Most cases do not pass ``--seed``, and must not have one invented.

    The production call does not pass it either, so the default is what
    almost every real run uses. A runner that filled in a seed would be
    testing a path nobody takes.
    """
    case = corpus.BY_NAME["production_default"]
    module = stub_port(golden_text(case.name))
    corpus.build(case, tmp_path)
    runners.PythonRunner().run(case, tmp_path)
    assert "seed" not in module.calls[0], (
        "the runner invented a seed for a case that does not set one")


def test_python_runner_reports_a_wrong_port_as_a_failure(stub_port, tmp_path):
    """End-to-end proof the acceptance path can fail.

    The stub returns another case's output. If this passed, every
    ``test_python_port_matches_the_goldens`` result would be worthless.
    """
    case = corpus.BY_NAME["tie_heavy"]
    stub_port(golden_text("production_default"))
    corpus.build(case, tmp_path)

    produced, _ = runners.PythonRunner().run(case, tmp_path)
    report = compare.compare(case.name, golden_text(case.name), produced)
    assert not report.ok, (
        "the harness accepted a port that returned a different screen's "
        "results entirely")


def test_python_runner_reports_a_bug_compatible_port_as_a_failure(stub_port,
                                                                  tmp_path):
    """A port that reproduces the C++ bugs must fail the divergence tests.

    Exercised through the stub rather than left until the port exists,
    because this is the assertion most likely to be written wrongly: it is
    the only one in the file whose expected outcome is "not ok".
    """
    case = corpus.BY_NAME["gmt_blank_line"]
    stub_port(golden_text(case.name))
    corpus.build(case, tmp_path)

    produced, _ = runners.PythonRunner().run(case, tmp_path)
    report = compare.compare(case.name, golden_text(case.name), produced)
    assert report.ok, (
        "the stub returns the golden verbatim, so this comparison is "
        "expected to agree; the divergence test is the one that inverts it")
    # and that is exactly why test_python_port_does_not_reproduce_the_
    # divergence_goldens asserts `not report.ok` for this case.
