"""Record goldens from the compiled RRA helper.

Run this while the C++ still exists::

    python tests/parity/record.py

It regenerates every case's input, runs the helper twice per case, and writes
the output as a golden only if both runs agree. The second run is not
redundant: it is the check that the oracle is deterministic, and therefore that
freezing its output is meaningful at all. RRA seeds its permutations with a
fixed ``PlantSeeds(123456)`` (``RRA.cpp:722,877``) so it should be -- but
``mageckGSEA`` seeded itself from the clock until recently, which is exactly
the failure this guard exists to catch.

``CONTRACT`` cases are skipped: they reach undefined behaviour in the C++, so
whatever it prints is an artifact of one build on one machine and must not
become a requirement.
"""

import argparse
import json
import pathlib
import platform
import subprocess
import sys
import time

# ``tests`` is deliberately not a package -- adding ``tests/__init__.py``
# would change how pytest imports every existing test module -- so this file
# has to work both as ``python tests/parity/record.py`` and as an import from
# the test suite, where pytest has already put ``tests/`` on sys.path.
if __package__ in (None, ""):
    sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
    from parity import compare, corpus, runners
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
    hits = corpus.build(case, directory)

    first, stdout = runner.run(case, directory)
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
    # The helper's own stdout is kept as evidence -- it records the per-stratum
    # permutation counts and any warnings -- but with trailing whitespace
    # stripped per line. RRA emits trailing spaces (for instance in the
    # "not found in the ranked list" warning), and committing them verbatim
    # would make `git diff --check` complain on every future touch of these
    # files. Nothing about the trailing spaces is information.
    (directory / "stdout.txt").write_text(
        "".join(line.rstrip() + "\n" for line in stdout.splitlines()))
    if not keep_inputs:
        (directory / "input.txt").unlink()

    rows, order = compare.parse(first)
    return {
        "why": case.why,
        "mode": case.mode.value,
        "alpha": case.alpha,
        "flags": list(case.flags),
        "uses_control": case.uses_control,
        "ngenes": case.ngenes,
        "n_rows": len(rows),
        "planted_hits": sorted(hits),
        "n_na_fdr": sum(1 for r in rows.values() if r["FDR"] is None),
        "min_p": min(r["p"] for r in rows.values()),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", default=str(GOLDENS), type=pathlib.Path,
                        help="directory to write goldens into")
    parser.add_argument("--binary", default="RRA",
                        help="name or path of the RRA helper")
    parser.add_argument("--case", action="append",
                        help="record only this case (repeatable)")
    args = parser.parse_args(argv)

    runner = runners.CppRunner(args.binary)
    if not runner.available():
        parser.error("%r is not on PATH. Goldens can only be recorded while "
                     "the compiled helper exists; run `pip install .` or "
                     "`make -C rra` first." % args.binary)

    selected = corpus.PARITY_CASES
    if args.case:
        unknown = set(args.case) - set(corpus.BY_NAME)
        if unknown:
            parser.error("unknown case(s): %s" % ", ".join(sorted(unknown)))
        selected = tuple(corpus.BY_NAME[n] for n in args.case)

    args.out.mkdir(parents=True, exist_ok=True)
    cases = {}
    for case in selected:
        if case.mode is corpus.Mode.CONTRACT:
            print("skip  %-24s (contract case: C++ behaviour is undefined)"
                  % case.name)
            continue
        started = time.time()
        cases[case.name] = record_case(case, runner, args.out)
        print("ok    %-24s %5.2fs  rows=%d"
              % (case.name, time.time() - started,
                 cases[case.name]["n_rows"]))

    manifest_path = args.out / MANIFEST
    existing = {}
    if manifest_path.exists():
        existing = json.loads(manifest_path.read_text()).get("cases", {})
    existing.update(cases)

    manifest_path.write_text(json.dumps({
        "recorded_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "recorded_by": "tests/parity/record.py",
        "oracle": {
            "implementation": "compiled RRA C++ helper",
            "fingerprint": runner.fingerprint(),
            "git_commit": _git("rev-parse", "HEAD"),
            "git_describe": _git("describe", "--tags", "--always", "--dirty"),
            # Scoped to the C++ trees on purpose. The question a future
            # reader has is "was this binary built from the recorded
            # commit?", and uncommitted changes under tests/ -- including
            # the harness itself -- do not affect that. A dirty rra/ does.
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
            "atol_fdr": compare.ATOL_FDR,
            "note": "lo_value and p are printed via %10.4e (five significant "
                    "figures) and FDR via %f, so the file itself bounds how "
                    "closely two implementations can be compared.",
        },
        "cases": existing,
    }, indent=2, sort_keys=True) + "\n")
    print("\nwrote %s (%d case(s))" % (manifest_path, len(existing)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
