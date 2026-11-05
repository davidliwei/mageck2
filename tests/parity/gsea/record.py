"""Record goldens from the compiled ``mageckGSEA`` helper.

Run this while the C++ still exists::

    python tests/parity/gsea/record.py

It regenerates every case's input, runs the helper twice per case, and writes
the output as a golden only if both runs agree. The second run is not
redundant: it is the check that the oracle is deterministic, and therefore that
freezing its output is meaningful at all. This helper in particular earned that
check -- it seeded its permutations from the clock until PR #52, so goldens
recorded from any earlier build would have been noise.

Malformed-input cases are recorded like any other. Their goldens are not a
requirement on the port -- the whole point of those cases is that the port must
differ -- but the C++ behaviour there is defined and reproducible, and
recording it is what keeps the premise of each fix checkable. ``CONTRACT``
cases, if any are added, are skipped: they reach undefined behaviour, so
whatever the helper prints is an artifact of one build on one machine.
"""

import argparse
import json
import pathlib
import platform
import subprocess
import sys
import time

# ``tests`` is deliberately not a package -- adding ``tests/__init__.py`` would
# change how pytest imports every existing test module -- so this file has to
# work both as ``python tests/parity/gsea/record.py`` and as an import from the
# test suite, where pytest has already put ``tests/`` on sys.path.
if __package__ in (None, ""):
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))
    from parity.gsea import compare, corpus, runners
else:
    from . import compare, corpus, runners

GOLDENS = pathlib.Path(__file__).parent / "goldens"
MANIFEST = "manifest.json"


def _git(*args):
    try:
        out = subprocess.run(("git",) + args, capture_output=True, text=True,
                             cwd=str(pathlib.Path(__file__).resolve().parent))
        return out.stdout.strip() if out.returncode == 0 else None
    except OSError:
        return None


def record_case(case, runner, root, keep_inputs=True):
    """Build, run twice, verify agreement, and write the golden.

    Returns a dict of metadata for the manifest.
    """
    directory = root / case.name
    directory.mkdir(parents=True, exist_ok=True)
    enriched = corpus.build(case, directory)

    first, log = runner.run(case, directory)
    second, _ = runner.run(case, directory)
    if first != second:
        report = compare.compare(case.name, first, second)
        raise RuntimeError(
            "the reference implementation is not deterministic on case %r, so "
            "its output cannot be used as a golden:\n%s"
            % (case.name, report.render()))

    # The runner writes its own scratch output into the case directory. The
    # golden is expected.txt, so the scratch copy is removed rather than
    # committed next to it as a confusing duplicate.
    _, scratch = runner.command(case, directory)
    if scratch.exists():
        scratch.unlink()

    (directory / "expected.txt").write_text(first)
    # The helper's own log is kept as evidence -- it records the pathway and
    # rank sizes it parsed, the seed it used, and any warnings -- but with
    # trailing whitespace stripped per line, since the progress dots and
    # warnings carry some and committing it would make `git diff --check`
    # complain on every future touch of these files.
    # Trailing blank lines go too. The helper writes an unconditional
    # `cerr << endl` after the pathway loop while only printing progress dots
    # when sorting, so an unsorted run's log ends in a bare newline, and
    # committing it trips `git diff --check` on every future touch.
    kept = [line.rstrip() for line in log.splitlines()]
    while kept and not kept[-1]:
        kept.pop()
    (directory / "stderr.txt").write_text(
        "".join(line + "\n" for line in kept))
    if not keep_inputs:
        (directory / "rank.txt").unlink()
        (directory / "sets.gmt").unlink()

    rows, order = compare.parse(first)
    return {
        "why": case.why,
        "mode": case.mode.value,
        "score_column": case.score_column,
        "perms": case.perms,
        "seed": case.seed,
        "flags": {
            "sort_by_p": case.sort_by_p,
            "skip_header": case.skip_header,
            "reverse": case.reverse,
            "pathway_name": case.pathway_name,
        },
        "ngenes": case.ngenes,
        "nsets": case.nsets,
        "n_rows": len(rows),
        "enriched_sets": sorted(enriched),
        "signal_expected": case.signal_expected,
        "n_degenerate": sum(1 for r in rows.values()
                            if r["ES"] == 0.0 and r["p_permutation"] == 1.0),
        "min_p_permutation": min(r["p_permutation"] for r in rows.values()),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", default=str(GOLDENS), type=pathlib.Path,
                        help="directory to write goldens into")
    parser.add_argument("--binary", default="mageckGSEA",
                        help="name or path of the mageckGSEA helper")
    parser.add_argument("--case", action="append",
                        help="record only this case (repeatable)")
    args = parser.parse_args(argv)

    runner = runners.CppRunner(args.binary)
    if not runner.available():
        parser.error("%r is not on PATH. Goldens can only be recorded while "
                     "the compiled helper exists; run `pip install .` or "
                     "`make -C gsea` first." % args.binary)

    selected = corpus.RECORDED_CASES
    if args.case:
        unknown = set(args.case) - set(corpus.BY_NAME)
        if unknown:
            parser.error("unknown case(s): %s" % ", ".join(sorted(unknown)))
        selected = tuple(corpus.BY_NAME[n] for n in args.case)

    args.out.mkdir(parents=True, exist_ok=True)
    cases = {}
    for case in selected:
        if case.mode is corpus.Mode.CONTRACT:
            print("skip  %-28s (contract case: C++ behaviour is undefined)"
                  % case.name)
            continue
        started = time.time()
        cases[case.name] = record_case(case, runner, args.out)
        print("ok    %-28s %5.2fs  rows=%-4d %s"
              % (case.name, time.time() - started,
                 cases[case.name]["n_rows"],
                 "(malformed input)" if case in corpus.MALFORMED_CASES else ""))

    manifest_path = args.out / MANIFEST
    existing = {}
    if manifest_path.exists():
        existing = json.loads(manifest_path.read_text()).get("cases", {})
    existing.update(cases)

    manifest_path.write_text(json.dumps({
        "recorded_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "recorded_by": "tests/parity/gsea/record.py",
        "oracle": {
            "implementation": "compiled mageckGSEA C++ helper",
            "fingerprint": runner.fingerprint(),
            "git_commit": _git("rev-parse", "HEAD"),
            "git_describe": _git("describe", "--tags", "--always", "--dirty"),
            # Scoped to the C++ trees on purpose. The question a future reader
            # has is "was this binary built from the recorded commit?", and
            # uncommitted changes under tests/ -- including the harness
            # itself -- do not affect that. A dirty gsea/ does.
            "oracle_source_dirty": bool(
                _git("status", "--porcelain", "--", "rra", "gsea")),
            "tree_dirty": bool(_git("status", "--porcelain")),
        },
        "platform": {
            "system": platform.system(),
            "machine": platform.machine(),
            "python": sys.version.split()[0],
        },
        "tolerances": {
            "rtol_printed": compare.RTOL_PRINTED,
            "atol_permutation": compare.ATOL_PERMUTATION,
            "atol_es": compare.ATOL_ES,
            "note": "every number is written through ostream << with default "
                    "flags, i.e. six significant digits, so the file itself "
                    "bounds how closely two implementations can be compared. "
                    "p_permutation is required to agree exactly: the "
                    "generator is seeded and its consumption order fully "
                    "determined, so there is no Monte Carlo noise to allow "
                    "for.",
            "es_is_platform_dependent": "The ES column is not reproducible to "
                    "the last bit, which is why atol_es exists. "
                    "getgseascore_fast computes decvalue*(ni-i) + sumincval, "
                    "which clang contracts into a fused multiply-add at -O3 "
                    "wherever FMA is in the baseline instruction set (it is "
                    "on the arm64 these were recorded on; it is not on "
                    "baseline x86-64). Rebuilding the same source with "
                    "-ffp-contract=off changes three of these goldens by "
                    "2.776e-17 and 5.551e-17, in gene sets whose score "
                    "cancels to zero. Read the platform block above before "
                    "concluding that a mismatch in that digit means anything.",
        },
        "cases": existing,
    }, indent=2, sort_keys=True) + "\n")
    print("\nwrote %s (%d case(s))" % (manifest_path, len(existing)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
