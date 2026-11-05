"""Rank-agreement statistics shared by both ports' comparators.

These live here rather than in either comparator because the question they
answer -- "would a user notice the difference?" -- is the same for RRA gene p
values and for GSEA pathway p values, and because an agreement statistic that
is subtly different between the two harnesses would make their reports
incomparable.

Nothing here pulls in SciPy. The suite already depends on it, but a tie-aware
Spearman and an inversion count are a dozen lines each, and writing them out
makes the tie handling -- which is the whole difficulty -- visible at the point
where it matters.
"""

import numpy as np


def average_ranks(x):
    """Ranks of ``x``, averaged within each run of equal values.

    Tie averaging is not a nicety here. Permutation p values are multiples of
    1/perms, so they pile up on a coarse grid and on the 0 floor; a tie-naive
    ranking would report a spurious disagreement between two implementations
    that put the same items on the floor in a different order.
    """
    x = np.asarray(x, dtype=float)
    order = np.argsort(x, kind="stable")
    ranks = np.empty(len(x), dtype=float)
    ranks[order] = np.arange(len(x), dtype=float)
    sorted_x = x[order]
    start = 0
    for i in range(1, len(x) + 1):
        if i == len(x) or sorted_x[i] != sorted_x[start]:
            ranks[order[start:i]] = ranks[order[start:i]].mean()
            start = i
    return ranks


def normalized_ranks(x):
    """Tie-averaged ranks scaled so 0 is the smallest value and 1 the largest.

    Used to express "these items rank near the top" in a way that does not
    depend on how many items there are.
    """
    x = np.asarray(x, dtype=float)
    if len(x) < 2:
        return np.zeros(len(x))
    return average_ranks(x) / (len(x) - 1)


def spearman(a, b):
    """Tie-aware Spearman correlation between two value vectors."""
    ra, rb = average_ranks(a), average_ranks(b)
    ra = ra - ra.mean()
    rb = rb - rb.mean()
    denom = np.sqrt((ra ** 2).sum() * (rb ** 2).sum())
    if denom == 0.0:
        return 1.0 if np.allclose(a, b) else 0.0
    return float((ra * rb).sum() / denom)


def inversions(a, b):
    """Count item pairs the two implementations order oppositely.

    Pairs that tie in either implementation are not counted: a tie carries no
    ordering claim, so breaking it differently is not an inversion. This is
    O(n^2) and the corpora are a few hundred items, which is the right trade
    for a diagnostic that has to be obviously correct.
    """
    n = len(a)
    count = 0
    for i in range(n):
        for j in range(i + 1, n):
            da, db = a[i] - a[j], b[i] - b[j]
            if da == 0.0 or db == 0.0:
                continue
            if (da > 0.0) != (db > 0.0):
                count += 1
    return count
