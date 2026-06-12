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


def test_package_imports():
    import mageck2  # noqa: F401
    from mageck2.version import __version__

    assert __version__


def test_rra_binary_on_path():
    assert shutil.which("RRA") is not None, "compiled RRA binary not found on PATH"


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
