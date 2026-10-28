"""Regression coverage for automatic UMI windows at read boundaries (#34)."""

import logging
from types import SimpleNamespace

import pytest

from mageck2.mageckCountIO import (
    mageckcount_determine_umi_location,
    mageckcount_search_variable_region,
)


@pytest.mark.parametrize(
    "sequences, expected",
    [
        ([base * 8 for base in "ACGT"] * 50, (0, 8)),
        (["AA" + base * 8 for base in "ACGT"] * 50, (2, 10)),
        ([base for base in "ACGT"] * 50, (0, 1)),
        (["AA" + base for base in "ACGT"] * 50, (2, 3)),
        ([base * 8 + "A" for base in "ACGT"] * 50, (0, 8)),
        ([base * 8 + "AA" for base in "ACGT"] * 50, (0, 8)),
        (["AAAA"] * 200, (-1, -1)),
        (["N"] * 200, (-1, -1)),
        ([""] * 200, (-1, -1)),
        ([], (-1, -1)),
        ([base * 8 for base in "ACGT"] * 49 + ["AAAAAAAAAA"] * 4, (0, 8)),
    ],
)
def test_variable_region_boundaries(sequences, expected):
    assert mageckcount_search_variable_region(sequences) == expected


def test_terminal_umi_is_selected_on_first_read(caplog):
    args = SimpleNamespace(umi="auto")
    with caplog.at_level(logging.INFO):
        mageckcount_determine_umi_location(
            args, None, {}, [base * 8 for base in "ACGT"] * 50
        )

    assert (args.umi, args.umi_start, args.umi_end) == ("firstpair", 0, 8)
    assert "UMI found in the first read. Position (after guide): 0-8" in caplog.text
