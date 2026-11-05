"""Differential-test harness for the C++ -> Python ports.

MAGeCK2 ships two compiled helpers, ``RRA`` and ``mageckGSEA``, and both are
being replaced by Python. They have been the reference implementations since
2013; their output is what thousands of published screens were analysed with,
so "the new code looks right" is not a sufficient standard. Each port has to be
shown to agree with the old implementation, case by case, on the paths users
actually exercise.

This package exists to record that agreement **while the C++ still exists**. A
``record.py`` runs the compiled helper over a fixed corpus and stores its
output as goldens; a ``compare.py`` then checks any candidate implementation
against those goldens. Once the C++ trees are deleted the goldens remain as the
oracle, and the last commit that contained the C++ remains in git history as a
way to regenerate them.

Four things make this harness trustworthy rather than decorative:

* It proves its own oracle. Each ``record.py`` runs the helper twice per case
  and refuses to write a golden if the runs disagree, so a non-deterministic
  reference cannot be silently frozen into a golden.
* It records provenance. Every golden set carries the git commit, the binary's
  SHA-256 and the build platform, so a golden can always be traced to the code
  that produced it.
* It separates parity from contract from known bugs. See :mod:`parity.modes`:
  some cases reach undefined behaviour in the C++, and some reach behaviour
  that is perfectly defined and wrong. Neither may be frozen into a
  requirement, and they are not the same thing.
* Its comparators are mutation-tested. A harness that cannot fail is not
  evidence, and "zero deviation everywhere" is also what a vacuous test
  reports, so each comparator has tests that inject a specific corruption and
  require it to be caught.

Layout
------

``modes`` and ``stats`` are shared. The RRA harness sits at the top level of
this package (``corpus``, ``compare``, ``runners``, ``record``, ``goldens``)
because it was here first; the GSEA harness sits in the :mod:`parity.gsea`
subpackage. Flattening that asymmetry would rename forty-odd committed golden
files for cosmetic gain, which is not worth doing while the ports are in
flight. Phase 4 deletes the C++ and is the natural moment to tidy it.
"""
