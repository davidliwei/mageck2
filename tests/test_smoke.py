"""Smoke tests for MAGeCK2.

These exercise the installed package end to end: the Python package imports,
the compiled RRA helper is on PATH, and `mageck2 test` runs the full ranking
pipeline on a small bundled count table and produces a well-formed
gene_summary file. They are intentionally lightweight so they can run on every
push in CI.
"""

import shutil
import subprocess
from pathlib import Path

DATA = Path(__file__).parent / "data" / "count_table.txt"
GMT = Path(__file__).parent / "data" / "pathways.gmt"


def test_package_imports():
    import mageck2  # noqa: F401
    from mageck2.version import __version__

    assert __version__


def test_rra_binary_on_path():
    assert shutil.which("RRA") is not None, "compiled RRA binary not found on PATH"


def test_mageckgsea_binary_on_path():
    assert (
        shutil.which("mageckGSEA") is not None
    ), "compiled mageckGSEA binary not found on PATH"


def test_mageck2_cli_runs():
    result = subprocess.run(
        ["mageck2", "--help"], capture_output=True, text=True
    )
    assert result.returncode == 0
    assert "count" in result.stdout and "test" in result.stdout


def test_mageck2_test_end_to_end(tmp_path):
    result = subprocess.run(
        [
            "mageck2", "test",
            "-k", str(DATA),
            "-t", "HL60.final,KBM7.final",
            "-c", "HL60.initial,KBM7.initial",
            "-n", "smoke",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr

    gene_summary = tmp_path / "smoke.gene_summary.txt"
    assert gene_summary.exists(), "gene_summary output not produced"

    lines = gene_summary.read_text().splitlines()
    header = lines[0].split("\t")
    assert "neg|score" in header
    assert "pos|score" in header
    assert len(lines) > 1, "no genes in gene_summary output"


def _run_test_step(tmp_path):
    """Produce a gene_summary the pathway step can consume."""
    result = subprocess.run(
        [
            "mageck2", "test",
            "-k", str(DATA),
            "-t", "HL60.final,KBM7.final",
            "-c", "HL60.initial,KBM7.initial",
            "-n", "smoke",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    gene_summary = tmp_path / "smoke.gene_summary.txt"
    assert gene_summary.exists(), "gene_summary output not produced"
    return gene_summary


def test_pathway_gsea_end_to_end(tmp_path):
    """`pathway --method gsea` (the default) runs via the mageckGSEA helper.

    Regression test for issue #10: the default GSEA method used to shell out to
    a mageckGSEA binary that was never built or shipped, so the run crashed on a
    missing intermediate file. The helper is now compiled and installed with
    MAGeCK2, so the full default path must produce a pathway summary.
    """
    gene_summary = _run_test_step(tmp_path)

    result = subprocess.run(
        [
            "mageck2", "pathway",
            "--gene-ranking", str(gene_summary),
            "--gmt-file", str(GMT),
            "-n", "pw",
            "--method", "gsea",
            "--permutation", "100",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr

    summary = tmp_path / "pw.pathway_summary.txt"
    assert summary.exists(), "pathway_summary output not produced"
    lines = summary.read_text().splitlines()
    assert len(lines) > 1, "no pathways in pathway_summary output"
    pathways = {line.split("\t")[0] for line in lines[1:]}
    assert {"PATHWAY_A", "PATHWAY_B", "PATHWAY_C"} <= pathways


def test_pathway_rra_end_to_end(tmp_path):
    """`pathway --method rra` runs via the RRA helper and writes a summary."""
    gene_summary = _run_test_step(tmp_path)

    result = subprocess.run(
        [
            "mageck2", "pathway",
            "--gene-ranking", str(gene_summary),
            "--gmt-file", str(GMT),
            "-n", "pw_rra",
            "--method", "rra",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "pw_rra.pathway_summary.txt").exists()


def test_gsea_skips_genes_absent_from_rank_file(tmp_path):
    """Regression: mageckGSEA must ignore pathway genes not in the rank file.

    They used to be inserted with a default index of 0 (via ``map::operator[]``)
    and treated as the top-ranked gene, silently inflating the enrichment score.
    Two pathways with the same real genes must therefore score identically even
    when one is padded with genes absent from the rank file, and a pathway made
    up entirely of absent genes must not crash (it has no overlap to score).
    """
    rank = tmp_path / "rank.txt"
    rank.write_text("id\tscore\n" + "".join(f"G{i}\t{i}\n" for i in range(1, 21)))

    gmt = tmp_path / "p.gmt"
    gmt.write_text(
        "REAL\tna\tG1\tG2\tG3\n"
        "PADDED\tna\tG1\tG2\tG3\tFAKE1\tFAKE2\tFAKE3\n"
        "ONLYFAKE\tna\tFAKE8\tFAKE9\n"
    )

    out = tmp_path / "out.txt"
    result = subprocess.run(
        [
            "mageckGSEA",
            "-c", "1",          # score column (0-indexed): id, score
            "-p", "100",
            "-g", str(gmt),
            "-r", str(rank),
            "-o", str(out),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr

    rows = {
        line.split("\t")[0]: line.split("\t")
        for line in out.read_text().splitlines()[1:]
    }
    # ES is column index 2. Padding with absent genes must not change it.
    assert rows["REAL"][2] == rows["PADDED"][2]
    # A pathway with no overlap scores 0 without crashing.
    assert float(rows["ONLYFAKE"][2]) == 0.0


def test_gsea_degenerate_pathway_is_not_maximally_significant(tmp_path):
    """Regression: a degenerate pathway must not report p-value 0.

    When a pathway covers the entire ranked list there is no enrichment to
    detect: ``getgseascore_fast`` returns a degenerate sentinel score of 0.
    Every permutation resamples the same full set and also scores 0, so the
    permutation test counted no exceedance and reported ``p_permutation=0`` --
    flagging the pathway as *maximally* significant, exactly backwards.
    ``getscoreandp`` must instead short-circuit and leave the p-values at 1.0.
    """
    # Four genes, no header line: every line is a gene, so the pathway below
    # covers the whole ranked list (n == pathway size).
    rank = tmp_path / "rank.txt"
    rank.write_text("A\t4\nB\t3\nC\t2\nD\t1\n")

    gmt = tmp_path / "p.gmt"
    gmt.write_text("ALL\tna\tA\tB\tC\tD\n")

    out = tmp_path / "out.txt"
    result = subprocess.run(
        [
            "mageckGSEA",
            "-c", "1",
            "-p", "100",
            "-g", str(gmt),
            "-r", str(rank),
            "-o", str(out),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr

    row = {
        line.split("\t")[0]: line.split("\t")
        for line in out.read_text().splitlines()[1:]
    }["ALL"]
    # ES is column 2; permutation p-value is column 4. A degenerate pathway
    # scores 0 and must NOT be reported as maximally significant (p=0).
    assert float(row[2]) == 0.0
    assert float(row[4]) == 1.0


def _gsea_es(tmp_path, rank_text, header_flag):
    """Run mageckGSEA on ``rank_text`` and return the enrichment score of PW."""
    rank = tmp_path / f"rank_{'H' if header_flag else 'plain'}.txt"
    rank.write_text(rank_text)
    gmt = tmp_path / "p.gmt"
    # G14..G18 hold the lowest scores, so they cluster at one end of the ranking
    # and produce a robustly non-zero enrichment score (~0.87) rather than a
    # near-zero value whose textual form differs across platforms.
    gmt.write_text("PW\tna\tG14\tG15\tG16\tG17\tG18\n")
    out = tmp_path / f"out_{'H' if header_flag else 'plain'}.txt"
    cmd = ["mageckGSEA", "-c", "1", "-p", "100",
           "-g", str(gmt), "-r", str(rank), "-o", str(out)]
    if header_flag:
        cmd.append("-H")
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    rows = {
        line.split("\t")[0]: line.split("\t")[2]
        for line in out.read_text().splitlines()[1:]
    }
    return rows["PW"]


def test_gsea_skip_header_matches_headerless(tmp_path):
    """Regression: ``-H`` must drop the rank-file header, not rank it as a gene.

    ``mageck2 pathway --method gsea`` passes a ``*.gene_summary.txt`` whose first
    row is a header. Without skipping it, mageckGSEA parsed the header as a
    phantom gene named ``id`` (score ``atof("neg|score") == 0``), inflating the
    gene count and shifting every pathway's ranking, percentile, and permutation
    math. A headered file read with ``-H`` must therefore score identically to
    the same data with no header line at all.
    """
    genes = "".join(f"G{i}\t{20 - i}\n" for i in range(1, 21))  # G1..G20, distinct scores
    headerless = _gsea_es(tmp_path, genes, header_flag=False)
    headered = _gsea_es(tmp_path, "id\tscore\n" + genes, header_flag=True)

    # -H makes the headered input score identically to the headerless data.
    assert headered == headerless
    # Sanity check that a real (non-degenerate) enrichment score was computed,
    # so the equality above is meaningful and not a pair of zeros.
    assert float(headered) > 0.5


def test_explicit_trim5_does_not_crash(tmp_path):
    """Regression test for mageck2-doc issue #1.

    An explicit numeric --trim-5 used to raise
    ``NameError: name 'candidate_trim5' is not defined`` before reading any
    reads. Here a 20 bp guide preceded by a 20 bp prefix is recovered with
    --trim-5 20, so the function must run to completion and return its
    (candidate_trim5_list, remaining_sequences) tuple.
    """
    from types import SimpleNamespace

    from mageck2.mageckCountIO import mageckcount_trim5_auto

    guide = "ACGTACGTACGTACGTACGT"          # 20 bp guide present in the library
    read = "TTTTTTTTTTTTTTTTTTTT" + guide    # 20 bp prefix to be trimmed, then the guide
    genedict = {guide: ("sg1", "geneA")}

    fastq = tmp_path / "tiny.fastq"
    fastq.write_text("@r1\n" + read + "\n+\n" + "I" * len(read) + "\n")

    args = SimpleNamespace(trim_5="20", count_n=False)
    candidate_trim5, remaining = mageckcount_trim5_auto(str(fastq), args, genedict)

    assert candidate_trim5 == [20]


def test_pg_pair_only_skips_header(tmp_path, caplog):
    import logging
    from types import SimpleNamespace

    from mageck2.mageckCount import mageckcount_checkargs, mageckcount_iter_pg_pair_file

    lib1 = tmp_path / "lib1.txt"
    lib1.write_text(
        "sgRNA\tsequence\tgene\n"
        "sg1\tACGTACGTACGTACGTACGT\tGENE1\n"
    )
    lib2 = tmp_path / "lib2.txt"
    lib2.write_text(
        "sgRNA\tsequence\tgene\n"
        "sg2\tTGCATGCATGCATGCATGCA\tGENE2\n"
    )
    pair_file = tmp_path / "guide_pair.txt"
    pair_file.write_text(
        "guide1_fix\tguide2_fix\n"
        "sg1\tsg2\n"
    )

    assert list(mageckcount_iter_pg_pair_file(str(pair_file))) == [
        (2, ["sg1", "sg2"])
    ]

    # tab-delimited ids that contain spaces must stay intact (not whitespace-split)
    spaced_pair = tmp_path / "guide_pair_spaced.txt"
    spaced_pair.write_text(
        "guide1_fix\tguide2_fix\n"
        "Non-Targeting Control 1\tNon-Targeting Control 2\n"
    )
    assert list(mageckcount_iter_pg_pair_file(str(spaced_pair))) == [
        (2, ["Non-Targeting Control 1", "Non-Targeting Control 2"])
    ]

    # .csv paired-guide files are split on commas, matching the count-table readers
    csv_pair = tmp_path / "guide_pair.csv"
    csv_pair.write_text(
        "guide1_fix,guide2_fix\n"
        "sg1,sg2\n"
    )
    assert list(mageckcount_iter_pg_pair_file(str(csv_pair))) == [
        (2, ["sg1", "sg2"])
    ]

    args = SimpleNamespace(
        fastq=["sample_R1.fastq"],
        fastq_2=["sample_R2.fastq"],
        list_seq=str(lib1),
        list_seq_2=str(lib2),
        reverse_complement=False,
        reverse_complement_2=False,
        sample_label="",
        pg_pair_only=str(pair_file),
        count_table=None,
        day0_label=None,
    )

    with caplog.at_level(logging.WARNING):
        mageckcount_checkargs(args)

    assert "guide1_fix" not in caplog.text
    assert "guide2_fix" not in caplog.text


def test_mle_bayes_reports_disabled(tmp_path, caplog):
    import logging
    from types import SimpleNamespace

    from mageck2.mlemageck import mageckmle_main

    args = SimpleNamespace(
        output_prefix=str(tmp_path / "bayes_check"),
        design_matrix="1,0;1,1",
        include_samples="sample1,sample2",
        beta_labels="baseline,treatment",
        bayes=True,
    )

    with caplog.at_level(logging.ERROR):
        try:
            mageckmle_main(parsedargs=args)
        except SystemExit as exc:
            assert exc.code == 1
        else:
            raise AssertionError("mle --bayes should exit with an error")

    assert "mle --bayes" in caplog.text
    assert "disabled" in caplog.text


def test_read_cnvdata_accepts_string_symbols(tmp_path):
    from mageck2.cnv_normalization import read_CNVdata

    cnv_file = tmp_path / "cnv.txt"
    cnv_file.write_text(
        "SYMBOL\tHL60\tKBM7\n"
        "GENE1\t0\t1\n"
        "GENE2\t2\t3\n"
    )

    cnv_arr, cell_dict, gene_dict = read_CNVdata(str(cnv_file), ["HL60"])

    assert gene_dict == {"GENE1": 0, "GENE2": 1}
    assert cell_dict == {"HL60": 0}
    assert cnv_arr.tolist() == [[1.0], [4.0]]


def _write_designmat(tmp_path, name, rows, header="Samples\tbaseline\tHL60\tKBM7"):
    p = tmp_path / name
    p.write_text(header + "\n" + "".join("\t".join(r) + "\n" for r in rows))
    return p


def _run_mle(tmp_path, designmat, prefix):
    return subprocess.run(
        [
            "mageck2", "mle",
            "-k", str(DATA),
            "-d", str(designmat),
            "-n", prefix,
            "--permutation-round", "1",
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )


# Rows for the bundled count table, in the order the model requires: the two
# baseline-only samples first, then one condition each.
_GOOD_ROWS = [
    ["HL60.initial", "1", "0", "0"],
    ["KBM7.initial", "1", "0", "0"],
    ["HL60.final", "1", "1", "0"],
    ["KBM7.final", "1", "0", "1"],
]


def test_mle_rejects_non_baseline_first_row(tmp_path):
    """Regression: the baseline sample must be the first row of the design matrix.

    ``DesignMatCache.save_record`` fits the model from ``design_mat[1:,1:]`` --
    row 0 is consumed as the reference and its own design row is never read. A
    matrix whose first row is a *condition* sample used to run to completion
    anyway: that sample was silently treated as a baseline and its condition
    discarded, biasing every beta measured against the reference. It must now be
    rejected before any counts are read.
    """
    rows = [_GOOD_ROWS[2]] + [_GOOD_ROWS[0], _GOOD_ROWS[1], _GOOD_ROWS[3]]
    designmat = _write_designmat(tmp_path, "badorder.txt", rows)

    result = _run_mle(tmp_path, designmat, "badorder")

    assert result.returncode != 0, "a non-baseline first row must not be accepted"
    combined = result.stdout + result.stderr
    assert "first row" in combined and "HL60.final" in combined
    assert not (tmp_path / "badorder.gene_summary.txt").exists()


def test_mle_rejects_designmat_without_baseline_sample(tmp_path):
    """Every sample carrying a condition leaves no reference to fit against."""
    rows = [
        ["HL60.initial", "1", "1", "0"],
        ["KBM7.initial", "1", "0", "1"],
        ["HL60.final", "1", "1", "0"],
        ["KBM7.final", "1", "0", "1"],
    ]
    designmat = _write_designmat(tmp_path, "nobase.txt", rows)

    result = _run_mle(tmp_path, designmat, "nobase")

    assert result.returncode != 0
    assert "no baseline sample" in result.stdout + result.stderr


def test_mle_rejects_non_unit_baseline_column(tmp_path):
    """The first column is consumed as the baseline term, so it must be all 1s."""
    rows = [
        ["HL60.initial", "1", "0", "0"],
        ["KBM7.initial", "0", "0", "0"],
        ["HL60.final", "1", "1", "0"],
        ["KBM7.final", "1", "0", "1"],
    ]
    designmat = _write_designmat(tmp_path, "badcol.txt", rows)

    result = _run_mle(tmp_path, designmat, "badcol")

    assert result.returncode != 0
    combined = result.stdout + result.stderr
    assert "baseline column" in combined and "KBM7.initial" in combined


def test_mle_accepts_valid_designmat_with_pooled_baselines(tmp_path):
    """The demo matrix -- two baseline-only samples pooled -- must still run.

    Guards the validation above against over-rejecting: multiple baseline rows
    are legal (they are pooled into one shared reference, see issue #23).
    """
    designmat = _write_designmat(tmp_path, "good.txt", _GOOD_ROWS)

    result = _run_mle(tmp_path, designmat, "good")

    assert result.returncode == 0, result.stderr
    gene_summary = tmp_path / "good.gene_summary.txt"
    assert gene_summary.exists()
    header = gene_summary.read_text().splitlines()[0].split("\t")
    assert "HL60|beta" in header and "KBM7|beta" in header


def test_documented_designmat_example_is_valid():
    """The example in `mageck2 mle --help` must satisfy the design-matrix rules.

    The help text advertised "1,1;1,0" as the simplest invocation, but that
    matrix leads with a condition sample, so validate_designmat rejects it --
    following the documentation verbatim produced an error. Extract whatever
    example the help text currently advertises and check it actually validates,
    so the two cannot drift apart again.
    """
    import re

    from mageck2.mledesignmat import parse_designmat, validate_designmat

    result = subprocess.run(
        ["mageck2", "mle", "--help"], capture_output=True, text=True
    )
    assert result.returncode == 0, result.stderr

    # collapse argparse's line wrapping before matching the quoted example
    helptext = " ".join(result.stdout.split())
    example = re.search(r'For example, "([0-9,;]+)"', helptext)
    assert example is not None, "no design-matrix example found in mle --help"

    desmat, _rows, _cols = parse_designmat(example.group(1))
    validate_designmat(desmat)  # raises SystemExit if the example is not usable


def test_duplicate_sequence_report_names_dropped_ids(tmp_path, caplog):
    """Library rows dropped for sharing a sequence must be identified.

    ``mageckcount_checklists`` keeps the first row for a given sequence and
    drops later ones. The report has to say which library it read, which IDs it
    dropped, and which ID now represents that sequence -- otherwise the user
    cannot connect the later "not in the first library" warning back to its
    cause. See issue #27.
    """
    import logging

    from mageck2.mageckCount import mageckcount_checklists

    lib = tmp_path / "lib1.txt"
    lib.write_text(
        "sgRNA\tsequence\tgene\n"
        "sg_keep\tACGTACGTACGTACGTACGT\tGENE1\n"
        "sg_dup\tACGTACGTACGTACGTACGT\tGENE1\n"
    )

    with caplog.at_level(logging.WARNING):
        genedict = mageckcount_checklists(str(lib), False)

    # the later row is dropped, as before
    assert "sg_keep" in genedict
    assert "sg_dup" not in genedict

    # ... but now the user is told which file, which ID, and its replacement
    assert "lib1.txt" in caplog.text
    assert "sg_dup" in caplog.text
    assert "sg_keep" in caplog.text


def test_clean_library_reports_no_duplicate_warning(tmp_path, caplog):
    """A library with no duplicated sequences must not warn.

    The report used to fire at WARNING level unconditionally, announcing
    "There are 0 sgRNAs with duplicated sequences." on every clean run, which
    trains users to ignore it. See issue #27.
    """
    import logging

    from mageck2.mageckCount import mageckcount_checklists

    lib = tmp_path / "clean.txt"
    lib.write_text(
        "sgRNA\tsequence\tgene\n"
        "sg1\tACGTACGTACGTACGTACGT\tGENE1\n"
        "sg2\tTGCATGCATGCATGCATGCA\tGENE2\n"
    )

    with caplog.at_level(logging.WARNING):
        genedict = mageckcount_checklists(str(lib), False)

    assert len(genedict) == 2
    assert caplog.text == ""


def test_pg_pair_only_warning_explains_duplicate_drop(tmp_path, caplog):
    """The --pg-pair-only warning must name the cause and the remedy.

    "sgRNA ID sg_dup not in the first library" is true but unactionable: the ID
    is missing because its sequence duplicates an earlier row, and the fix is to
    use the surviving representative ID. See issue #27.
    """
    import logging
    from types import SimpleNamespace

    from mageck2.mageckCount import mageckcount_checkargs

    lib1 = tmp_path / "lib1.txt"
    lib1.write_text(
        "sgRNA\tsequence\tgene\n"
        "sg_keep\tACGTACGTACGTACGTACGT\tGENE1\n"
        "sg_dup\tACGTACGTACGTACGTACGT\tGENE1\n"
    )
    lib2 = tmp_path / "lib2.txt"
    lib2.write_text(
        "sgRNA\tsequence\tgene\n"
        "sg2\tTGCATGCATGCATGCATGCA\tGENE2\n"
    )
    pair_file = tmp_path / "pairs.txt"
    pair_file.write_text(
        "guide1\tguide2\n"
        "sg_dup\tsg2\n"
    )

    args = SimpleNamespace(
        fastq=["sample_R1.fastq"],
        fastq_2=["sample_R2.fastq"],
        list_seq=str(lib1),
        list_seq_2=str(lib2),
        reverse_complement=False,
        reverse_complement_2=False,
        sample_label="",
        pg_pair_only=str(pair_file),
        count_table=None,
        day0_label=None,
    )

    with caplog.at_level(logging.WARNING):
        mageckcount_checkargs(args)

    # Isolate the pair-file warning from the library-loading one, which also
    # mentions both IDs.
    pair_messages = [r.getMessage() for r in caplog.records if "pairs.txt" in r.getMessage()]
    assert len(pair_messages) == 1
    message = pair_messages[0]

    assert "sg_dup" in message
    assert "sg_keep" in message      # the representative to use instead
    assert "duplicate" in message.lower()


def _revcomp(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _write_paired_guide_fixture(tmp_path, unmatched_first_reads=0, unmatched_second_reads=0,
                                r2_prefix="", lib2_revcomp=False):
    """A two-construct dual-guide library, one construct allowed by the pair file.

    r2_prefix pads read 2 ahead of the second guide, so the guide no longer sits
    at offset 0 and the extraction window has to be told where it is.
    lib2_revcomp stores the second library reverse-complemented relative to the
    reads, which is what --reverse-complement-2 exists to undo.
    """
    guide1 = "ACGTACGTACGTACGTACGT"
    guide1b = "AACCGGTTAACCGGTTAACC"  # in the library, never sequenced
    guide2a = "TTGACCTAGGCATTGACCAA"
    guide2b = "GGGGCCCCAAAATTTTGGGG"

    # a guide that is its own reverse complement matches in both orientations, so
    # it cannot tell a correct --reverse-complement-2 from a missing one
    assert _revcomp(guide2a) != guide2a
    assert _revcomp(guide2b) != guide2b

    (tmp_path / "lib1.txt").write_text(
        "sgRNA\tsequence\tgene\n"
        f"sg1\t{guide1}\tGENE1\n"
        f"sg1b\t{guide1b}\tGENE1\n"
    )
    lib2 = (guide2a, guide2b)
    if lib2_revcomp:
        lib2 = tuple(_revcomp(g) for g in lib2)
    (tmp_path / "lib2.txt").write_text(
        "sgRNA\tsequence\tgene\n"
        f"sg2a\t{lib2[0]}\tGENE2\n"
        f"sg2b\t{lib2[1]}\tGENE3\n"
    )
    # only sg1+sg2a is a real construct; sg1+sg2b is a recombinant to be filtered
    (tmp_path / "pairs.txt").write_text(
        "guide1\tguide2\n"
        "sg1\tsg2a\n"
    )

    def fastq(reads):
        return "".join(
            f"@r{i}\n{seq}\n+\n{'I' * len(seq)}\n" for i, seq in enumerate(reads)
        )

    r1 = [guide1] * 7
    r2 = [r2_prefix + guide2a] * 4 + [r2_prefix + guide2b] * 3
    # read pairs whose first guide matches nothing in the library
    r1 += ["TTTTTTTTTTTTTTTTTTTT"] * unmatched_first_reads
    r2 += [r2_prefix + guide2a] * unmatched_first_reads
    # read pairs whose second guide matches nothing in the second library
    r1 += [guide1] * unmatched_second_reads
    r2 += [r2_prefix + "CCCCCCCCCCCCCCCCCCCC"] * unmatched_second_reads
    (tmp_path / "r1.fastq").write_text(fastq(r1))
    (tmp_path / "r2.fastq").write_text(fastq(r2))


def _run_paired_guide_count(tmp_path, extra_args=(), start2="0", end2="20", pair_only=True):
    return subprocess.run(
        [
            "mageck2", "count",
            "-l", "lib1.txt",
            "--list-seq-2", "lib2.txt",
            "--pairguide", "secondpair",
            "--pg-start-2", start2,
            "--pg-end-2", end2,
            *(("--pg-pair-only", "pairs.txt") if pair_only else ()),
            "--trim-5", "0",
            "--sample-label", "s1",
            "--fastq", "r1.fastq",
            "--fastq-2", "r2.fastq",
            "-n", "pg",
            *extra_args,
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
    )


def test_paired_guide_count_end_to_end(tmp_path):
    """Baseline: the allowed pair is counted, the filtered pair is not."""
    _write_paired_guide_fixture(tmp_path)
    result = _run_paired_guide_count(tmp_path)
    assert result.returncode == 0, result.stderr

    rows = (tmp_path / "pg.pg_count.txt").read_text().splitlines()
    assert rows[0].split("\t")[:2] == ["sgRNA1_sgRNA2", "Gene1_Gene2"]

    counts = {r.split("\t")[0]: r.split("\t")[2] for r in rows[1:]}
    assert counts == {"sg1_sg2a": "4"}


def test_pg_unmapped_file_records_filtered_pairs(tmp_path):
    """Pairs excluded from pg_count.txt must be recorded, not silently dropped.

    mageckcount_printpgdict is currently called with ounmappedfile=None, so a
    pair rejected by --pg-pair-only (or whose second guide never matched) leaves
    no trace anywhere, and the user cannot tell a filtered pair from one that
    was never sequenced. See issue #27.
    """
    _write_paired_guide_fixture(tmp_path)
    result = _run_paired_guide_count(tmp_path, extra_args=["--unmapped-to-file"])
    assert result.returncode == 0, result.stderr

    pg_unmapped = tmp_path / "pg.pg_unmapped.txt"
    assert pg_unmapped.exists(), "no pair-level unmapped file was written"

    rows = pg_unmapped.read_text().splitlines()
    filtered = {r.split("\t")[0]: r.split("\t")[2] for r in rows[1:]}
    assert filtered == {"sg1_sg2b": "3"}


def test_pg_unmapped_file_tolerates_unmatched_first_reads(tmp_path):
    """Read pairs whose first guide never matched must not break the writer.

    In secondpair mode the second-read window is recorded even when the first
    read matched nothing, under a None key. That key was previously invisible
    because the pair-level unmapped file was never written; once it is, the row
    label None + '_' + seq raises TypeError. Those pairs identify no guide 1 and
    are already accounted for in the read-level unmapped file, so they belong in
    neither pair-level table. See issue #27.
    """
    _write_paired_guide_fixture(tmp_path, unmatched_first_reads=2)
    result = _run_paired_guide_count(tmp_path, extra_args=["--unmapped-to-file"])
    assert result.returncode == 0, result.stderr

    rows = (tmp_path / "pg.pg_unmapped.txt").read_text().splitlines()
    filtered = {r.split("\t")[0]: r.split("\t")[2] for r in rows[1:]}
    assert filtered == {"sg1_sg2b": "3"}

    # the unmatched first reads are still reported, at the read level
    read_unmapped = (tmp_path / "pg.unmapped.txt").read_text()
    assert "TTTTTTTTTTTTTTTTTTTT" in read_unmapped


def test_pg_unmapped_rows_match_their_header(tmp_path):
    """Every row of pg_unmapped.txt must mean what the header says.

    A pair whose second guide matched no entry in --list-seq-2 has no second ID
    and no second gene, but the row was serialized as sequences under an
    "sgRNA1_sgRNA2 / Gene1_Gene2" header: column 1 held two raw sequences and
    column 2 held an sgRNA ID joined to a sequence. Anything parsing this file
    reads a 20-mer where an ID belongs. Report the known first guide by ID and
    gene, keep the unmatched second sequence as its only available identifier,
    and mark the missing gene explicitly.
    """
    _write_paired_guide_fixture(tmp_path, unmatched_second_reads=2)
    result = _run_paired_guide_count(tmp_path, extra_args=["--unmapped-to-file"])
    assert result.returncode == 0, result.stderr

    rows = [r.split("\t") for r in (tmp_path / "pg.pg_unmapped.txt").read_text().splitlines()]
    assert rows[0][:2] == ["sgRNA1_sgRNA2", "Gene1_Gene2"]

    unmapped = {r[0]: (r[1], r[2]) for r in rows[1:]}

    # excluded by --pg-pair-only: both guides are known, so both are named
    assert unmapped["sg1_sg2b"] == ("GENE1_GENE3", "3")

    # second guide unknown: first guide by ID and gene, second by sequence
    assert unmapped["sg1_CCCCCCCCCCCCCCCCCCCC"] == ("GENE1_None", "2")


def _write_duplicate_sequence_pg_fixture(tmp_path):
    """A dual-guide library indexed by construct, so one position-1 ID is a
    duplicate-sequence casualty -- the shape reported in issue #15.

    sg1 and sg1_dup carry the same sequence, so sg1_dup is dropped on load and
    sg1 represents that sequence. The pair file, written against the original
    construct nomenclature, names sg1_dup for the second construct.
    """
    guide1 = "ACGTACGTACGTACGTACGT"
    guide1b = "AACCGGTTAACCGGTTAACC"  # in the library, never sequenced
    guide2a = "TGCATGCATGCATGCATGCA"
    guide2b = "GGGGCCCCAAAATTTTGGGG"

    (tmp_path / "lib1.txt").write_text(
        "sgRNA\tsequence\tgene\n"
        f"sg1\t{guide1}\tGENE1\n"
        f"sg1_dup\t{guide1}\tGENE1\n"
        f"sg1b\t{guide1b}\tGENE1\n"
    )
    (tmp_path / "lib2.txt").write_text(
        "sgRNA\tsequence\tgene\n"
        f"sg2a\t{guide2a}\tGENE2\n"
        f"sg2b\t{guide2b}\tGENE3\n"
    )
    (tmp_path / "pairs.txt").write_text(
        "guide1\tguide2\n"
        "sg1\tsg2a\n"
        "sg1_dup\tsg2b\n"
    )

    def fastq(reads):
        return "".join(
            f"@r{i}\n{seq}\n+\n{'I' * len(seq)}\n" for i, seq in enumerate(reads)
        )

    (tmp_path / "r1.fastq").write_text(fastq([guide1] * 7))
    (tmp_path / "r2.fastq").write_text(fastq([guide2a] * 4 + [guide2b] * 3))


def test_pg_pair_only_matches_pairs_naming_dropped_ids(tmp_path):
    """A pair file naming a dropped duplicate ID must still match its construct.

    --pg-pair-only is matched by ID, but counting resolves a sequence to the one
    surviving ID, so a pair naming a dropped ID could never match and the whole
    construct was excluded from pg_count.txt -- despite its sequence pair having
    been counted correctly. Users had to rewrite the pair file by hand. The
    remapping is mechanical, so do it. See issue #28.
    """
    _write_duplicate_sequence_pg_fixture(tmp_path)
    result = _run_paired_guide_count(tmp_path)
    assert result.returncode == 0, result.stderr

    rows = (tmp_path / "pg.pg_count.txt").read_text().splitlines()
    counts = {r.split("\t")[0]: r.split("\t")[2] for r in rows[1:]}

    # sg1_dup was dropped on load, so its construct is labelled with sg1
    assert counts == {"sg1_sg2a": "4", "sg1_sg2b": "3"}


def test_pg_pair_only_warning_reports_resolution_not_exclusion(tmp_path, caplog):
    """The startup warning must describe what now happens, not what used to.

    It told the user the pair would be excluded and to edit the pair file. Both
    are obsolete once the ID is resolved automatically. See issue #28.
    """
    import logging
    from types import SimpleNamespace

    from mageck2.mageckCount import mageckcount_checkargs

    _write_duplicate_sequence_pg_fixture(tmp_path)

    args = SimpleNamespace(
        fastq=["r1.fastq"],
        fastq_2=["r2.fastq"],
        list_seq=str(tmp_path / "lib1.txt"),
        list_seq_2=str(tmp_path / "lib2.txt"),
        reverse_complement=False,
        reverse_complement_2=False,
        sample_label="",
        pg_pair_only=str(tmp_path / "pairs.txt"),
        count_table=None,
        day0_label=None,
    )

    with caplog.at_level(logging.WARNING):
        mageckcount_checkargs(args)

    pair_messages = [r.getMessage() for r in caplog.records if "pairs.txt" in r.getMessage()]
    assert len(pair_messages) == 1
    message = pair_messages[0]

    assert "sg1_dup" in message
    assert "sg1" in message
    assert "resolved" in message.lower()
    assert "excluded" not in message.lower()


def test_pg_pair_only_logs_alias_summary(tmp_path):
    """Resolution must never be silent.

    On a construct-indexed library most of the pair file is remapped, and the
    user has to be able to see how many entries were affected and what they
    became. See issue #28.
    """
    _write_duplicate_sequence_pg_fixture(tmp_path)
    result = _run_paired_guide_count(tmp_path)
    assert result.returncode == 0, result.stderr

    assert "1 of 2 records" in result.stderr
    assert "sg1_dup -> sg1" in result.stderr


def test_pg_pair_only_reports_colliding_resolved_pairs(tmp_path):
    """Two entries that resolve to the same pair must be reported, not merged
    in silence.

    Uniqueness has to hold at the sequence-pair level, not the ID-pair level.
    When both libraries contain duplicated sequences, two distinct pair-file
    entries can resolve onto one row; their reads are then indistinguishable and
    are summed. That is unavoidable -- but the user has to be told. See #28.
    """
    guide1 = "ACGTACGTACGTACGTACGT"
    guide1b = "AACCGGTTAACCGGTTAACC"
    guide2a = "TGCATGCATGCATGCATGCA"

    (tmp_path / "lib1.txt").write_text(
        "sgRNA\tsequence\tgene\n"
        f"sg1\t{guide1}\tGENE1\n"
        f"sg1_dup\t{guide1}\tGENE1\n"
        f"sg1b\t{guide1b}\tGENE1\n"
    )
    (tmp_path / "lib2.txt").write_text(
        "sgRNA\tsequence\tgene\n"
        f"sg2a\t{guide2a}\tGENE2\n"
        f"sg2a_dup\t{guide2a}\tGENE2\n"
    )
    (tmp_path / "pairs.txt").write_text(
        "guide1\tguide2\n"
        "sg1\tsg2a\n"
        "sg1_dup\tsg2a_dup\n"
    )
    reads = [(guide1, guide2a)] * 5
    (tmp_path / "r1.fastq").write_text(
        "".join(f"@r{i}\n{a}\n+\n{'I' * len(a)}\n" for i, (a, _) in enumerate(reads))
    )
    (tmp_path / "r2.fastq").write_text(
        "".join(f"@r{i}\n{b}\n+\n{'I' * len(b)}\n" for i, (_, b) in enumerate(reads))
    )

    result = _run_paired_guide_count(tmp_path)
    assert result.returncode == 0, result.stderr

    # both entries name the same sequence pair, so there is one row carrying all reads
    rows = (tmp_path / "pg.pg_count.txt").read_text().splitlines()
    assert {r.split("\t")[0]: r.split("\t")[2] for r in rows[1:]} == {"sg1_sg2a": "5"}

    assert "resolve to the same pair" in result.stderr
    assert "sg1_dup" in result.stderr and "sg2a_dup" in result.stderr


def test_pg_pair_only_does_not_remap_ids_that_survive_loading(tmp_path):
    """An ID re-admitted under a different sequence must not be aliased away.

    The library loader guards duplicate IDs against what it has already loaded,
    so an ID dropped for duplicating an earlier sequence can be accepted later
    under a sequence of its own. It then holds a stale entry in the dropped map,
    and resolving pair-file IDs through that map unconditionally would redirect a
    perfectly valid ID to an unrelated representative -- allowing the wrong
    sequence pair and filtering the intended one. Found by review of PR #31.
    """
    shared = "ACGTACGTACGTACGTACGT"   # claimed by sgA first
    own = "TTTTAAAACCCCGGGGTTTT"      # sgX's own sequence, admitted later
    other = "AACCGGTTAACCGGTTAACC"
    guide2 = "TGCATGCATGCATGCATGCA"

    (tmp_path / "lib1.txt").write_text(
        "sgRNA\tsequence\tgene\n"
        f"sgA\t{shared}\tGENE1\n"
        f"sgX\t{shared}\tGENE1\n"   # dropped here: sequence already used by sgA
        f"sgX\t{own}\tGENE2\n"      # but re-admitted here, with a sequence of its own
        f"sgB\t{other}\tGENE3\n"
    )
    (tmp_path / "lib2.txt").write_text(
        "sgRNA\tsequence\tgene\n"
        f"sg2\t{guide2}\tGENE4\n"
    )
    (tmp_path / "pairs.txt").write_text(
        "guide1\tguide2\n"
        "sgX\tsg2\n"
    )
    reads = 4
    (tmp_path / "r1.fastq").write_text(
        "".join(f"@r{i}\n{own}\n+\n{'I' * len(own)}\n" for i in range(reads))
    )
    (tmp_path / "r2.fastq").write_text(
        "".join(f"@r{i}\n{guide2}\n+\n{'I' * len(guide2)}\n" for i in range(reads))
    )

    result = _run_paired_guide_count(tmp_path)
    assert result.returncode == 0, result.stderr

    rows = (tmp_path / "pg.pg_count.txt").read_text().splitlines()
    counts = {r.split("\t")[0]: r.split("\t")[2] for r in rows[1:]}

    # sgX is a real, loaded sgRNA: its pair must be allowed under its own name
    assert counts == {"sgX_sg2": "4"}


def test_ambiguous_duplicate_id_is_left_unresolved(tmp_path, caplog):
    """An ID that names two different duplicated sequences cannot be aliased.

    Each duplicate-sequence row overwrote the ID's entry in the dropped map, so
    the ID resolved to whichever representative happened to come last -- making
    which sequence pair passes --pg-pair-only depend on library row order. Such
    an ID is genuinely ambiguous: leave it unresolved and say so. Found by
    review of PR #31.
    """
    import logging

    from mageck2.mageckCount import mageckcount_checklists

    seq_s = "ACGTACGTACGTACGTACGT"
    seq_t = "TGCATGCATGCATGCATGCA"
    header = "sgRNA\tsequence\tgene\n"
    claimed = f"sgA\t{seq_s}\tGENE1\nsgB\t{seq_t}\tGENE2\n"

    # sgX duplicates seq_s and seq_t, each already represented by an earlier row
    one_order = tmp_path / "lib_sx_st.txt"
    one_order.write_text(header + claimed + f"sgX\t{seq_s}\tGENE1\nsgX\t{seq_t}\tGENE2\n")
    other_order = tmp_path / "lib_st_sx.txt"
    other_order.write_text(header + claimed + f"sgX\t{seq_t}\tGENE2\nsgX\t{seq_s}\tGENE1\n")

    with caplog.at_level(logging.WARNING):
        for lib in (one_order, other_order):
            dropped = {}
            genedict = mageckcount_checklists(str(lib), False, dropped=dropped)

            assert "sgX" not in genedict     # neither row was loaded
            assert "sgX" not in dropped      # and it resolves to neither sgA nor sgB

    ambiguity = [r.getMessage() for r in caplog.records if "sgX" in r.getMessage()]
    assert len(ambiguity) == 2               # reported for each library, not silent
    assert "sgA" in ambiguity[0] and "sgB" in ambiguity[0]


def _write_umi_fixture(tmp_path, umi_on="second", trailing=True, nread=200):
    """A single-guide library plus reads carrying a random 8bp UMI.

    umi_on selects which read holds the UMI; trailing=False makes read 1 end
    exactly at the guide, so nothing follows it to search.
    """
    import random

    rng = random.Random(11)
    const = "CTTGTGGAAAGGACGAAACA"
    guides = [
        ("sg%d" % i, "".join(rng.choice("ACGT") for _ in range(20)), "GENE%d" % (i % 3))
        for i in range(1, 9)
    ]
    (tmp_path / "lib1.txt").write_text(
        "sgRNA\tsequence\tgene\n" + "".join(f"{a}\t{b}\t{c}\n" for a, b, c in guides)
    )
    r1, r2 = [], []
    for i in range(nread):
        guide = guides[i % len(guides)][1]
        umi = "".join(rng.choice("ACGT") for _ in range(8))
        if umi_on == "first":
            r1.append(guide + umi)
            # variable read 2: the search finds a start here but never an end,
            # which is the (0, -1) failure value the guard has to reject
            r2.append("".join(rng.choice("ACGT") for _ in range(20)))
        else:
            r1.append(guide + (const if trailing else ""))
            r2.append(umi + const)

    def fastq(reads):
        return "".join(f"@r{i}\n{s}\n+\n{'I' * len(s)}\n" for i, s in enumerate(reads))

    (tmp_path / "r1.fastq").write_text(fastq(r1))
    (tmp_path / "r2.fastq").write_text(fastq(r2))


def _run_count(tmp_path, extra):
    return subprocess.run(
        ["mageck2", "count", "-l", "lib1.txt", "--trim-5", "0", "--sample-label", "s1",
         "--fastq", "r1.fastq", "--fastq-2", "r2.fastq", "-n", "out", *extra],
        cwd=tmp_path, capture_output=True, text=True,
    )


def test_pairguide_auto_is_rejected(tmp_path):
    """--pairguide auto never located a second guide, so it is no longer offered.

    The variable-region search it used finds a column run whose base composition
    looks random. That describes a UMI, not a guide drawn from a fixed library,
    so the window it returned was truncated or empty and pg_count.txt came out
    header-only with a zero exit code. Reject the choice outright rather than
    keep a mode that cannot succeed. See issue #29.
    """
    _write_paired_guide_fixture(tmp_path)
    result = _run_paired_guide_count(tmp_path, extra_args=["--pairguide", "auto"])

    assert result.returncode != 0
    assert "auto" in result.stderr
    assert "invalid choice" in result.stderr


def test_pairguide_secondpair_requires_its_window(tmp_path):
    """Omitting --pg-start-2/--pg-end-2 must stop the run, once and up front.

    The window was checked inside the per-read loop, which logged an error for
    every read -- one line per read on a real FASTQ -- then carried on to slice
    an empty string and exit 0 with a header-only pg_count.txt. See issue #29.
    """
    _write_paired_guide_fixture(tmp_path)
    result = subprocess.run(
        ["mageck2", "count", "-l", "lib1.txt", "--list-seq-2", "lib2.txt",
         "--pairguide", "secondpair", "--trim-5", "0", "--sample-label", "s1",
         "--fastq", "r1.fastq", "--fastq-2", "r2.fastq", "-n", "pg"],
        cwd=tmp_path, capture_output=True, text=True,
    )

    assert result.returncode != 0
    assert "--pg-start-2" in result.stderr and "--pg-end-2" in result.stderr
    # reported once, not once per read
    assert len([l for l in result.stderr.splitlines() if "ERROR" in l]) == 1
    assert not (tmp_path / "pg.pg_count.txt").exists()


def test_pairguide_firstpair_requires_its_window(tmp_path):
    """Same guard for the first-read layout, which names different options."""
    _write_paired_guide_fixture(tmp_path)
    result = subprocess.run(
        ["mageck2", "count", "-l", "lib1.txt", "--list-seq-2", "lib2.txt",
         "--pairguide", "firstpair", "--trim-5", "0", "--sample-label", "s1",
         "--fastq", "r1.fastq", "--fastq-2", "r2.fastq", "-n", "pg"],
        cwd=tmp_path, capture_output=True, text=True,
    )

    assert result.returncode != 0
    assert "--pg-start" in result.stderr and "--pg-end" in result.stderr
    assert len([l for l in result.stderr.splitlines() if "ERROR" in l]) == 1


def test_umi_secondpair_requires_its_window(tmp_path):
    """The same missing-window flood reaches --umi, which shares the code path."""
    _write_umi_fixture(tmp_path)
    result = _run_count(tmp_path, ["--umi", "secondpair"])

    assert result.returncode != 0
    assert "--umi-start-2" in result.stderr and "--umi-end-2" in result.stderr
    assert len([l for l in result.stderr.splitlines() if "ERROR" in l]) == 1


def test_umi_auto_fails_loudly_when_no_variable_region_is_found(tmp_path):
    """A failed search must not be accepted as a result.

    The second-read search was accepted when *either* bound was set, so the
    failure value (0, -1) passed as success. Reads were then sliced [0:-1] --
    the whole read but its last base -- and every read yielded a distinct
    "UMI", filling umi_count.txt with rows that look entirely plausible. Only
    a search that set both bounds has found anything. See issue #29.
    """
    _write_umi_fixture(tmp_path, umi_on="first")
    result = _run_count(tmp_path, ["--umi", "auto"])

    assert result.returncode != 0, "a failed UMI search was reported as success"
    assert "ailed" in result.stderr
    assert not (tmp_path / "out.umi_count.txt").exists()


def test_umi_auto_reports_reads_with_nothing_after_the_guide(tmp_path):
    """Reads that end at the guide leave no sequence to search.

    The column-frequency table is then empty and max() over its keys raised an
    uncaught ValueError, so the run died with a traceback rather than a message
    naming the problem. See issue #29.
    """
    _write_umi_fixture(tmp_path, trailing=False)
    result = _run_count(tmp_path, ["--umi", "auto"])

    assert "Traceback" not in result.stderr
    assert "ValueError" not in result.stderr


def test_unmatched_reads_do_not_accumulate_a_none_key(tmp_path):
    """Read pairs whose first guide matched nothing must not be retained.

    In secondpair mode the second-read window was recorded for every read,
    including those where read 1 matched no guide, under a None key. Nothing
    ever reads that entry -- mageckcount_printpgdict skips it -- but it grows
    with one entry per distinct unmatched window, so the memory held is sized by
    the *unmapped* fraction of the FASTQ rather than by the library. See
    issue #29.
    """
    import random

    from mageck2.mageckCountIO import mageckcount_processonefile

    rng = random.Random(3)

    def rnd(n):
        return "".join(rng.choice("ACGT") for _ in range(n))

    guides = [("sg%d" % i, rnd(20), "GENE%d" % i) for i in range(1, 6)]
    lib = tmp_path / "lib1.txt"
    lib.write_text("".join(f"{a}\t{b}\t{c}\n" for a, b, c in guides))

    matched, unmatched = 100, 2000
    r1 = [guides[i % len(guides)][1] for i in range(matched)]
    r1 += [rnd(20) for _ in range(unmatched)]   # in the library's alphabet, but absent from it
    r2 = [rnd(20) for _ in range(matched + unmatched)]

    def fastq(reads):
        return "".join(f"@r{i}\n{s}\n+\n{'I' * len(s)}\n" for i, s in enumerate(reads))

    (tmp_path / "r1.fastq").write_text(fastq(r1))
    (tmp_path / "r2.fastq").write_text(fastq(r2))

    class _Args:
        trim_5 = "0"
        count_n = False
        sgrna_len = 0
        unmapped_to_file = False
        test_run = False
        reverse_complement = False
        umi = "secondpair"
        umi_start = -1
        umi_end = -1
        umi_start_2 = 0
        umi_end_2 = 20
        pairguide = "secondpair"
        pg_start = -1
        pg_end = -1
        pg_start_2 = 0
        pg_end_2 = 20

    genedict = {seq: (sgid, gene) for sgid, seq, gene in guides}
    ctab, ctab_umi = {}, {}
    mageckcount_processonefile(
        str(tmp_path / "r1.fastq"), _Args(), ctab, ctab_umi, genedict, {},
        str(tmp_path / "r2.fastq"), adjust=True,
    )

    assert None not in ctab_umi
    # what is kept is bounded by the library, not by the unmapped read count
    assert set(ctab_umi) <= set(genedict)
    assert sum(len(v) for v in ctab_umi.values()) <= matched


def _pg_counts(tmp_path):
    rows = (tmp_path / "pg.pg_count.txt").read_text().splitlines()
    assert rows[0].split("\t")[:2] == ["sgRNA1_sgRNA2", "Gene1_Gene2"]
    return {r.split("\t")[0]: r.split("\t")[2] for r in rows[1:]}


def test_wrong_second_read_window_yields_no_pairs(tmp_path):
    """A misplaced --pg-start-2/--pg-end-2 is the quiet failure behind #15.

    count.txt fills normally from read 1 while pg_count.txt comes out with only
    a header, because the window slices sequence that matches nothing in
    --list-seq-2. Nothing about the exit code says so, so pin the shape of the
    failure: empty under the wrong window, correct under the right one, same
    reads either way. See issue #29.
    """
    _write_paired_guide_fixture(tmp_path, r2_prefix="CCCCC")

    missed = _run_paired_guide_count(tmp_path, start2="0", end2="20")
    assert missed.returncode == 0, missed.stderr
    assert _pg_counts(tmp_path) == {}

    # read 1 matched throughout: the single-guide table is unaffected
    single = (tmp_path / "pg.count.txt").read_text().splitlines()
    assert [r.split("\t")[2] for r in single[1:] if r.startswith("sg1\t")] == ["7"]

    found = _run_paired_guide_count(tmp_path, start2="5", end2="25")
    assert found.returncode == 0, found.stderr
    assert _pg_counts(tmp_path) == {"sg1_sg2a": "4"}


def test_reverse_complement_2_matches_only_the_orientation_it_is_given(tmp_path):
    """--reverse-complement-2 flips the library, not the read.

    It must be required when the library is stored opposite the reads, and must
    break the match when it is not -- a flag that appeared to work either way
    would mean the guides were self-complementary, not that the option works.
    See issue #29.
    """
    forward = tmp_path / "forward"
    forward.mkdir()
    _write_paired_guide_fixture(forward, lib2_revcomp=False)
    assert _run_paired_guide_count(forward).returncode == 0
    assert _pg_counts(forward) == {"sg1_sg2a": "4"}

    assert _run_paired_guide_count(forward, extra_args=["--reverse-complement-2"]).returncode == 0
    assert _pg_counts(forward) == {}, "the flag must not match a library already in read orientation"

    flipped = tmp_path / "flipped"
    flipped.mkdir()
    _write_paired_guide_fixture(flipped, lib2_revcomp=True)
    assert _run_paired_guide_count(flipped).returncode == 0
    assert _pg_counts(flipped) == {}, "a reverse-complemented library must not match without the flag"

    assert _run_paired_guide_count(flipped, extra_args=["--reverse-complement-2"]).returncode == 0
    assert _pg_counts(flipped) == {"sg1_sg2a": "4"}


def test_pg_min_read_excludes_pairs_below_the_threshold(tmp_path):
    """--pg-min-read drops sparse pairs; the boundary value is kept.

    The threshold is applied to the total across samples, and the comparison is
    >=, so a pair with exactly --pg-min-read reads must survive. See issue #29.
    """
    _write_paired_guide_fixture(tmp_path)

    assert _run_paired_guide_count(tmp_path, extra_args=["--pg-min-read", "4"]).returncode == 0
    assert _pg_counts(tmp_path) == {"sg1_sg2a": "4"}, "a pair at the threshold must be kept"

    assert _run_paired_guide_count(tmp_path, extra_args=["--pg-min-read", "5"]).returncode == 0
    assert _pg_counts(tmp_path) == {}

    # the default is 3, which the 4-read pair clears on its own
    assert _run_paired_guide_count(tmp_path).returncode == 0
    assert _pg_counts(tmp_path) == {"sg1_sg2a": "4"}


def test_all_pairs_are_reported_when_pg_pair_only_is_omitted(tmp_path):
    """Without --pg-pair-only every observed combination is reported.

    The filtered run keeps only sg1+sg2a; without the option the recombinant
    sg1+sg2b must appear too, so the option is doing the filtering rather than
    some other part of the path silently dropping the pair. See issue #29.
    """
    _write_paired_guide_fixture(tmp_path)

    assert _run_paired_guide_count(tmp_path, pair_only=False).returncode == 0
    assert _pg_counts(tmp_path) == {"sg1_sg2a": "4", "sg1_sg2b": "3"}

    assert _run_paired_guide_count(tmp_path, pair_only=True).returncode == 0
    assert _pg_counts(tmp_path) == {"sg1_sg2a": "4"}
