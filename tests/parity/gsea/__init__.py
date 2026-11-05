"""Differential-test harness for the ``mageckGSEA`` C++ -> Python port.

``mageckGSEA`` is smaller than ``RRA`` -- 577 lines of logic over the same
vendored math and RNG -- but it is not simply easier. Three things about it
make the harness shape different:

* **Its scoring depends on how the rank file was parsed.** ``parserankfile``
  has two unrelated branches. With ``-c 0`` the score is derived from a gene's
  line number; with ``-c N`` the genes are re-sorted by the chosen column and
  the score is derived from the new position. The two branches do not even use
  the same offset -- one counts lines from 1, the other indexes from 0 -- and
  both divide by a line count that is one larger than the number of genes. A
  port that "cleans this up" changes every score in the file.

* **The enrichment score is stored as a 32-bit float.** ``gs_ret::gscore`` is
  a ``float`` while everything computing it is a ``double``, and the
  permutation test compares at that truncated precision
  (``calculategsea.cpp:431``). The truncation is therefore not a display
  detail: it decides permutation p values. It also turns out to be protective
  -- see the next point.

* **The enrichment score is not reproducible to the last bit, and the C++
  does not agree with itself about it.** Each step is
  ``decvalue*(ni-i) + sumincval``, which clang at ``-O3`` contracts into a
  fused multiply-add wherever FMA is in the baseline instruction set. Rounding
  once instead of twice changes the last bit, which is invisible after the
  float32 truncation *except* where the true score cancels to zero: such a set
  scores 0 on a non-contracted build and 5.55112e-17 on a contracted one.
  Measured by rebuilding this source with ``-ffp-contract=off``, which moves
  three of the twenty goldens. So the port is not asked to match those digits;
  ``compare.ATOL_ES`` says what it is asked for instead, and the permutation p
  values -- which the truncation protects -- are still required to be exact.

* **Output row order is genuinely unspecified in places.** Pathways are sorted
  with ``std::sort`` on ``(p_permutation, -ES)``, which is unstable, and with a
  coarse permutation grid whole blocks of pathways tie on both keys. The
  ``Ranking`` column is that row position. The FDR values are not affected --
  the step-down pass collapses each tied run to one value -- but ``Ranking``
  and row order are, so the comparator cannot require them.

The modules mirror the RRA harness: ``corpus`` builds inputs, ``runners``
invokes an implementation, ``compare`` diffs two outputs column by column, and
``record`` freezes the C++ output as goldens.
"""
