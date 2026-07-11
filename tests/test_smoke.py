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
