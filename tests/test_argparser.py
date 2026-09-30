"""Tests for the argument parser's construction.

``build_parser()`` exists so the parser can be built without parsing anything:
the documentation build imports it to generate the command-line reference, and
an import that consumed ``sys.argv`` or exited would break that build. These
tests pin that property down, since nothing else in the suite would notice if
parser construction regained a side effect.
"""

import argparse

import pytest

from mageck2.argsParser import build_parser

#: Every subcommand the CLI is expected to expose.
SUBCOMMANDS = ["count", "test", "pathway", "plot", "mle"]


def test_build_parser_returns_parser_without_reading_argv(monkeypatch):
    """Building the parser must not look at the process arguments."""
    monkeypatch.setattr("sys.argv", ["mageck2", "--nonsense-that-would-fail"])

    parser = build_parser()

    assert isinstance(parser, argparse.ArgumentParser)


def test_build_parser_exposes_every_subcommand():
    parser = build_parser()

    actions = [a for a in parser._actions if isinstance(a, argparse._SubParsersAction)]
    assert len(actions) == 1, "expected exactly one subparser group"

    assert sorted(actions[0].choices) == sorted(SUBCOMMANDS)


def test_parsing_a_count_command_still_works():
    """The split must leave parsing behaviour untouched."""
    args = build_parser().parse_args(
        ["count", "-l", "lib.txt", "-n", "out", "--fastq", "a.fastq"]
    )

    assert args.subcmd == "count"
    assert args.list_seq == "lib.txt"
    assert args.output_prefix == "out"
    assert args.fastq == ["a.fastq"]


def test_pairguide_choices_exclude_auto():
    """`--pairguide auto` was removed in 0.3.0; see davidliwei/mageck2#32."""
    parser = build_parser()
    count = [
        a for a in parser._actions if isinstance(a, argparse._SubParsersAction)
    ][0].choices["count"]

    pairguide = [a for a in count._actions if a.dest == "pairguide"][0]
    assert list(pairguide.choices) == ["none", "firstpair", "secondpair"]


def test_building_the_parser_twice_is_safe():
    """The docs build may construct the parser more than once per process."""
    assert build_parser() is not build_parser()


def test_help_does_not_escape_as_an_exception(capsys):
    """`--help` exits cleanly rather than raising something unexpected."""
    with pytest.raises(SystemExit) as excinfo:
        build_parser().parse_args(["--help"])

    assert excinfo.value.code == 0
    assert "sgRNA" in capsys.readouterr().out
