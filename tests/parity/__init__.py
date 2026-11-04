"""Differential-test harness for the RRA C++ -> Python port.

The port replaces a compiled helper that has been the reference implementation
since 2013. Its output is what thousands of published screens were analysed
with, so "the new code looks right" is not a sufficient standard: the port has
to be shown to agree with the old implementation, case by case, on the paths
users actually exercise.

This package exists to record that agreement **while the C++ still exists**.
``record.py`` runs the compiled helper over a fixed corpus and stores its
output as goldens; ``compare.py`` then checks any candidate implementation
against those goldens. Once the C++ trees are deleted the goldens remain as the
oracle, and the last commit that contained the C++ remains in git history as a
way to regenerate them.

Three things make this harness trustworthy rather than decorative:

* It proves its own oracle. ``record.py`` runs the helper twice per case and
  refuses to write a golden if the runs disagree, so a non-deterministic
  reference cannot be silently frozen into a golden.
* It records provenance. Every golden set carries the git commit, the binary's
  SHA-256 and the build platform, so a golden can always be traced to the code
  that produced it.
* It separates parity from contract. Some cases reach genuinely undefined
  behaviour in the C++ (see ``Mode.CONTRACT`` in ``corpus.py``). Freezing a
  golden there would enshrine whatever the compiler happened to do, so those
  cases assert a stated contract for the new implementation instead.
"""
