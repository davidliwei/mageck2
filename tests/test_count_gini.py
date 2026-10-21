"""Small targeted libraries should still produce count summaries."""

import csv
import subprocess

import pytest

from mageck2.mageckCountIO import mageckcount_gini


@pytest.mark.parametrize("counts", [[0], [1], [100], [0, 0], [0, 0, 0], []])
def test_degenerate_gini(counts):
    assert mageckcount_gini(counts) == 0.0


@pytest.mark.parametrize("counts, expected", [([1, 1], 0.0), ([0, 10], 1.0), ([1, 2, 3], 1 / 3)])
def test_multiple_guide_gini(counts, expected):
    assert mageckcount_gini(counts) == pytest.approx(expected)


def test_count_single_guide_fastq(tmp_path):
    guide = "ACGTACGTACGTACGTACGT"
    library = tmp_path / "library.tsv"
    library.write_text("sgRNA\tsequence\tgene\nsg1\t" + guide + "\tGENE1\n")
    reads = tmp_path / "reads.fastq"
    reads.write_text("".join(f"@read{i}\n{guide}\n+\n{'I' * len(guide)}\n" for i in range(10)))

    result = subprocess.run(
        ["mageck2", "count", "-l", str(library), "--fastq", str(reads),
         "--trim-5", "0", "--sample-label", "sample", "-n", str(tmp_path / "out")],
        cwd=tmp_path, capture_output=True, text=True,
    )

    assert result.returncode == 0, result.stderr
    with (tmp_path / "out.count.txt").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 1
    assert rows[0]["sgRNA"] == "sg1"
    assert int(rows[0]["sample"]) == 10
    with (tmp_path / "out.countsummary.txt").open() as handle:
        summary = list(csv.DictReader(handle, delimiter="\t"))
    assert len(summary) == 1
    assert float(summary[0]["GiniIndex"]) == 0.0
