"""Tests for the Python GSEA implementation's own contract.

``tests/test_gsea_parity.py`` holds it to agreement with the compiled helper,
column by column, against goldens recorded from that helper. This file covers
what parity cannot: the environment switch, the array API, the seed bounds,
and the input-handling bugs the port is required *not* to reproduce.

The last of those is the reason the two files are separate. A golden recorded
from the C++ is the right oracle for behaviour the port must copy and the
wrong one for behaviour it must fix, and the fixes are asserted here as
properties -- a malformed file must give the same answer as its repair -- not
as expected numbers.
"""

import os
import shutil
import subprocess
from pathlib import Path

import pytest

from mageck2 import gsea, pathwayFunc

DATA = Path(__file__).parent / "data" / "count_table.txt"
GMT = Path(__file__).parent / "data" / "pathways.gmt"


# --------------------------------------------------------------------------
# The switch
# --------------------------------------------------------------------------

def test_the_compiled_helper_is_still_the_default(monkeypatch):
    monkeypatch.delenv(pathwayFunc.GSEA_IMPL_ENV, raising=False)
    assert pathwayFunc.gsea_implementation() == "cpp"


@pytest.mark.parametrize("value", ["python", "Python", " PYTHON "])
def test_the_switch_selects_the_port(monkeypatch, value):
    monkeypatch.setenv(pathwayFunc.GSEA_IMPL_ENV, value)
    assert pathwayFunc.gsea_implementation() == "python"


def test_an_unknown_implementation_is_refused(monkeypatch):
    """A typo must stop the run rather than silently pick one.

    Defaulting to either implementation on an unrecognised value would mean a
    misspelled switch produced results from the one the user did not ask for,
    with nothing in the log to say so.
    """
    monkeypatch.setenv(pathwayFunc.GSEA_IMPL_ENV, "pyton")
    with pytest.raises(SystemExit) as excinfo:
        pathwayFunc.gsea_implementation()
    assert excinfo.value.code != 0


# --------------------------------------------------------------------------
# The seed
# --------------------------------------------------------------------------

def test_the_default_seed_is_the_helpers():
    """Both implementations must start from the same place by default.

    If these diverged, every permutation p value would differ between the
    backends on identical input, and the parity suite would be comparing two
    different runs rather than two implementations.
    """
    assert gsea.DEFAULT_SEED == 123456
    assert gsea.MAX_SEED == 2147483646


@pytest.mark.parametrize("seed", [0, -1, 2147483647, 2147483648])
def test_an_out_of_range_seed_is_refused(seed, tmp_path):
    """The helper rejects these rather than folding them.

    A seed of 0 or of MODULUS reduces to 0, where the generator is a fixed
    point: every draw is 0 and every permutation is the same one. Folding
    silently would produce confident p values from a degenerate null.
    """
    gmt, rank = _tiny_input(tmp_path)
    with pytest.raises(ValueError) as excinfo:
        gsea.enrich_from_file(str(gmt), str(rank), seed=seed, perms=5)
    assert "seed" in str(excinfo.value)


# --------------------------------------------------------------------------
# The array API
# --------------------------------------------------------------------------

def _ranked(n=200):
    """A ranked list with a planted block of co-enriched genes."""
    genes = ["G%04d" % i for i in range(n)]
    scores = [float(i) for i in range(n)]
    return genes, scores


def test_the_api_takes_sequences_and_returns_an_object():
    """The point of the port: no temporary files, no parsing the output back."""
    genes, scores = _ranked()
    pathways = {"TOP": genes[:20], "BOTTOM": genes[-20:],
                "SPREAD": genes[::10]}

    result = gsea.enrich(pathways, genes, scores, perms=100)

    assert len(result) == 3
    assert set(result.names) == {"TOP", "BOTTOM", "SPREAD"}
    rows = list(result)
    assert all(set(r) >= {"pathway", "es", "p", "p_permutation", "fdr"}
               for r in rows)
    # to_text() must stay parseable as the historical format, since that is
    # what keeps the recorded goldens usable after the C++ is gone.
    text = result.to_text()
    assert text.startswith(gsea.GSEAResult.HEADER)
    assert len(text.strip().splitlines()) == 4


def test_only_the_order_of_the_scores_matters():
    """The walk weights come from rank, not from score magnitude.

    ``parserankfile`` sorts by ``(value, name)`` and then overwrites every
    value with one derived from its position, so any strictly monotone
    transform of the scores is the same input. Worth pinning: a caller
    passing log fold changes where another passes p values must get the same
    answer, and a future optimization that used the magnitudes would be a
    silent change in meaning.
    """
    genes, scores = _ranked()
    pathways = {"TOP": genes[:20], "SPREAD": genes[::10]}

    plain = gsea.enrich(pathways, genes, scores, perms=50)
    # Strictly increasing, wildly different scale and offset.
    stretched = [s * 1e6 + 7.5 for s in scores]
    same = gsea.enrich(pathways, genes, stretched, perms=50)

    assert plain.to_text() == same.to_text()


def test_the_api_is_reproducible_without_being_asked():
    genes, scores = _ranked()
    pathways = {"TOP": genes[:20], "SPREAD": genes[::10]}
    first = gsea.enrich(pathways, genes, scores, perms=100)
    second = gsea.enrich(pathways, genes, scores, perms=100)
    assert first.to_text() == second.to_text()


def test_a_different_seed_moves_the_permutation_p_but_not_the_score():
    """The enrichment score is deterministic; only its p value is sampled."""
    genes, scores = _ranked()
    pathways = {"TOP": genes[:20], "SPREAD": genes[::10]}
    a = gsea.enrich(pathways, genes, scores, perms=200, seed=123456)
    b = gsea.enrich(pathways, genes, scores, perms=200, seed=999983)

    assert list(a.names) == list(b.names)
    assert list(a.es) == list(b.es)
    assert list(a.p) == list(b.p)
    assert list(a.p_permutation) != list(b.p_permutation)


def test_mismatched_input_lengths_are_refused():
    genes, scores = _ranked()
    with pytest.raises(ValueError):
        gsea.enrich({"S": genes[:5]}, genes, scores[:-1], perms=5)


def test_a_set_whose_genes_are_all_absent_is_reported_not_crashed():
    """Real GMT files contain sets with no overlap at all.

    Regression boundary for #46: absent genes used to be inserted at position
    0 by ``map::operator[]`` and treated as the top-ranked gene, which turned
    a set with no overlap into the most enriched one in the file.
    """
    genes, scores = _ranked()
    pathways = {"REAL": genes[:10], "ABSENT": ["NOSUCHGENE%d" % i
                                               for i in range(10)]}
    result = gsea.enrich(pathways, genes, scores, perms=50)

    rows = {r["pathway"]: r for r in result}
    assert rows["ABSENT"]["es"] == 0.0
    assert rows["ABSENT"]["p"] == 1.0
    assert rows["REAL"]["es"] > 0.0


# --------------------------------------------------------------------------
# The input-handling bugs the port must not reproduce
# --------------------------------------------------------------------------

def _tiny_input(tmp_path, description="generated_gene_set"):
    gmt = tmp_path / "sets.gmt"
    rank = tmp_path / "rank.txt"
    gmt.write_text("SET\t%s\tG0001\tG0002\tG0003\n" % description)
    # Three columns, because -c 2 counts columns after the gene name and so
    # needs a third field to read -- as the helper does, which errors on a
    # narrower file rather than defaulting the score.
    rank.write_text("id\tnum\tscore\n" + "".join(
        "G%04d\t5\t%.3f\n" % (i, i * 0.1) for i in range(1, 21)))
    return gmt, rank


def test_a_description_with_spaces_does_not_become_gene_names(tmp_path):
    """GMT is tab-delimited; the helper tokenizes on whitespace.

    The second column is a free-text description, so its second and later
    words are read as members by the helper. That inflates ``Size`` on
    essentially every real GMT file -- MSigDB descriptions are prose -- and
    when one of those words happens to match a gene symbol it is scored, which
    moves ``ES`` and ``p``. The bundled ``pathways.gmt`` shows it: its sets
    hold five genes and the helper reports ``Size`` 7.
    """
    one_word = tmp_path / "one"
    many_words = tmp_path / "many"
    one_word.mkdir()
    many_words.mkdir()
    gmt_a, rank_a = _tiny_input(one_word, description="setdescription")
    gmt_b, rank_b = _tiny_input(many_words,
                                description="set description with words")

    a = gsea.enrich_from_file(str(gmt_a), str(rank_a), perms=20)
    b = gsea.enrich_from_file(str(gmt_b), str(rank_b), perms=20)

    assert list(a.sizes) == [3], "the set has three genes"
    assert a.to_text() == b.to_text(), (
        "the description changed the result, so its words are being read as "
        "members")


def test_trailing_whitespace_does_not_add_a_phantom_member(tmp_path):
    clean = tmp_path / "clean"
    dirty = tmp_path / "dirty"
    clean.mkdir()
    dirty.mkdir()
    gmt_a, rank_a = _tiny_input(clean)
    gmt_b, rank_b = _tiny_input(dirty)
    gmt_b.write_text(gmt_b.read_text().replace("\n", " \t \n"))

    a = gsea.enrich_from_file(str(gmt_a), str(rank_a), perms=20)
    b = gsea.enrich_from_file(str(gmt_b), str(rank_b), perms=20)
    assert a.to_text() == b.to_text()


def test_a_blank_gmt_line_is_not_a_gene_set(tmp_path):
    """The helper's worst input bug, and the quietest.

    A blank line becomes a set named with the empty string, which then matches
    the default ``--pathway-name`` of ``""``. The helper takes its
    single-set branch, reports that one empty set, and silently skips every
    real set in the file -- exiting 0, with a plausible-looking output file.
    """
    clean = tmp_path / "clean"
    blank = tmp_path / "blank"
    clean.mkdir()
    blank.mkdir()
    gmt_a, rank_a = _tiny_input(clean)
    gmt_b, rank_b = _tiny_input(blank)
    gmt_b.write_text("\n" + gmt_b.read_text() + "\n\n")

    a = gsea.enrich_from_file(str(gmt_a), str(rank_a), perms=20)
    b = gsea.enrich_from_file(str(gmt_b), str(rank_b), perms=20)

    assert "" not in set(b.names), "a blank line produced a gene set"
    assert a.to_text() == b.to_text()


@pytest.mark.parametrize("which", ["sets.gmt", "rank.txt"])
def test_a_missing_final_newline_changes_nothing(tmp_path, which):
    """The helper drops the last line of either file.

    In the rank file that is worse than losing a row: the dropped gene also
    lowers the line count every weight is divided by, so one missing newline
    perturbs the whole ranking.
    """
    whole = tmp_path / "whole"
    cut = tmp_path / "cut"
    whole.mkdir()
    cut.mkdir()
    gmt_a, rank_a = _tiny_input(whole)
    gmt_b, rank_b = _tiny_input(cut)
    target = cut / which
    target.write_text(target.read_text().rstrip("\n"))

    a = gsea.enrich_from_file(str(gmt_a), str(rank_a), perms=20)
    b = gsea.enrich_from_file(str(gmt_b), str(rank_b), perms=20)
    assert a.to_text() == b.to_text()


# --------------------------------------------------------------------------
# End to end
# --------------------------------------------------------------------------

def _gene_summary(tmp_path):
    proc = subprocess.run(
        ["mageck2", "test", "-k", str(DATA),
         "-t", "HL60.final,KBM7.final", "-c", "HL60.initial,KBM7.initial",
         "-n", "out"],
        cwd=str(tmp_path), capture_output=True, text=True)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    return tmp_path / "out.gene_summary.txt"


@pytest.mark.skipif(
    shutil.which("mageck2") is None or shutil.which("mageckGSEA") is None,
    reason="needs an installed mageck2 and the compiled helper")
def test_both_implementations_give_the_same_pathway_summary(tmp_path):
    """A real ``mageck2 pathway`` run, both ways.

    Every column must agree exactly, ``Size`` included. It was once exempt:
    the bundled GMT's descriptions are three words each, and the helper used
    to tokenize the GMT on whitespace and count two of those words as set
    members, so it reported 7 where the port reported 5. That was a bug in the
    helper, fixed in the C++ rather than imitated in the port, and the two now
    agree on this fixture byte for byte.
    """
    gene_summary = _gene_summary(tmp_path)

    summaries = {}
    for impl in ("cpp", "python"):
        workdir = tmp_path / impl
        workdir.mkdir()
        env = dict(os.environ, MAGECK2_GSEA_IMPL=impl)
        proc = subprocess.run(
            ["mageck2", "pathway", "--gene-ranking", str(gene_summary),
             "--gmt-file", str(GMT), "-n", "pw", "--method", "gsea",
             "--permutation", "500"],
            cwd=str(workdir), env=env, capture_output=True, text=True)
        assert proc.returncode == 0, proc.stdout + proc.stderr
        if impl == "python":
            assert "Python GSEA implementation" in proc.stdout + proc.stderr, (
                "the switch did not take effect, so this compared the "
                "compiled helper against itself")
        summaries[impl] = (workdir / "pw.pathway_summary.txt").read_text()

    cpp_rows = [ln.split("\t") for ln in summaries["cpp"].splitlines()]
    py_rows = [ln.split("\t") for ln in summaries["python"].splitlines()]
    assert len(cpp_rows) == len(py_rows)

    for cpp_row, py_row in zip(cpp_rows[1:], py_rows[1:]):
        assert cpp_row[0] == py_row[0], "row order differs"
        assert len(cpp_row) == len(py_row), (
            "%s has %d columns one way and %d the other"
            % (cpp_row[0], len(cpp_row), len(py_row)))
        for column, (a, b) in enumerate(zip(cpp_row, py_row)):
            assert a == b, (
                "column %d (%s) of %s differs: %r vs %r"
                % (column, cpp_rows[0][column], cpp_row[0], a, b))

    # Agreement alone would also hold if both sides were wrong in the same
    # way, so check the sizes against the GMT itself. Every description here
    # is multi-word, which is exactly what a whitespace-tokenizing parser
    # turns into two extra members per set -- so this pins the column that
    # used to be exempt to ground truth, not just to the other backend.
    members = {}
    for line in GMT.read_text().splitlines():
        fields = line.split("\t")
        if len(fields) <= 2:
            continue
        members[fields[0]] = len(fields) - 2
        assert " " in fields[1], (
            "%s no longer has a multi-word description, so this test can no "
            "longer tell a whitespace-tokenizing parser from a correct one"
            % fields[0])
    assert members, "the bundled GMT has no usable rows"

    num = cpp_rows[0].index("num")
    for cpp_row in cpp_rows[1:]:
        assert int(cpp_row[num]) == members[cpp_row[0]], (
            "%s: reported size %s, but the GMT lists %d members"
            % (cpp_row[0], cpp_row[num], members[cpp_row[0]]))
