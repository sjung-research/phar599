"""Analysis 4 - pixel-intensity histograms with identical, integer-aligned bins.

Bins are built over the full representable range ``[0, 2**bits - 1]``. Each bin is a
half-open interval of ``width`` consecutive integer values, so no bin contains more
discrete levels than another (non-integer bin widths on integer data produce
artificial comb patterns). If the last bin extends past the maximum value it holds
fewer valid levels; ``n_levels_in_bin`` records this and ``density_per_level``
corrects for it.

Counting uses ``np.bincount`` on the exact integer values, so totals are exact.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .stats import select


def integer_bin_edges(max_value: int, n_bins: int) -> np.ndarray:
    """Integer edges from 0 with equal integer width covering ``[0, max_value]``.

    ``width = ceil((max_value + 1) / n_bins)``; only as many bins as needed to reach
    ``max_value`` are returned, so the actual count can be slightly below ``n_bins``
    (e.g. 4096 levels / 300 requested -> width 14, 293 bins).
    """
    if n_bins < 1 or max_value < 0:
        raise ValueError("n_bins must be >= 1 and max_value >= 0")
    n_levels = int(max_value) + 1
    width = -(-n_levels // int(n_bins))  # ceil
    n_eff = -(-n_levels // width)
    return np.arange(0, width * n_eff + 1, width, dtype=np.int64)


def level_counts(values: np.ndarray, max_value: int, mask: np.ndarray | None = None) -> np.ndarray:
    """Exact count of every integer level 0..max_value."""
    v = select(values, mask)
    if v.size and not np.issubdtype(v.dtype, np.integer):
        raise TypeError("histograms require integer raw data")
    if v.size and (v.min() < 0 or v.max() > max_value):
        raise ValueError(f"values outside [0, {max_value}]")
    return np.bincount(v.astype(np.int64, copy=False), minlength=int(max_value) + 1)


def binned_counts(levels: np.ndarray, edges: np.ndarray) -> np.ndarray:
    padded = np.zeros(int(edges[-1]), dtype=np.int64)
    padded[: levels.size] = levels
    return np.add.reduceat(padded, edges[:-1].astype(np.intp))


def histogram_table(levels: np.ndarray, edges: np.ndarray, max_value: int) -> pd.DataFrame:
    counts = binned_counts(levels, edges)
    total = int(levels.sum())
    left, right = edges[:-1], edges[1:]
    n_levels = np.minimum(right, max_value + 1) - left
    frac = counts / total if total else np.full(counts.shape, np.nan)
    return pd.DataFrame(
        {
            "bin_left": left,
            "bin_right_exclusive": right,
            "n_levels_in_bin": n_levels,
            "count": counts,
            "fraction": frac,
            "density_per_level": frac / n_levels,
        }
    )


def survival_table(levels: np.ndarray, thresholds: np.ndarray) -> pd.DataFrame:
    """Fraction of pixels with value >= T for each integer threshold T."""
    total = int(levels.sum())
    tail = np.cumsum(levels[::-1])[::-1]  # tail[t] = count(v >= t)
    t = np.asarray(thresholds, dtype=np.int64)
    t = t[(t >= 0) & (t < levels.size)]
    return pd.DataFrame({"threshold": t, "fraction_at_or_above": tail[t] / total if total else np.nan})


def quantiles_from_levels(levels: np.ndarray, qs: tuple[float, ...] = (0.5, 0.9, 0.99, 0.998)) -> dict:
    """Lower (inverted-CDF) quantiles from exact level counts, for histogram annotation."""
    cdf = np.cumsum(levels) / max(int(levels.sum()), 1)
    return {f"q{q:g}": int(np.searchsorted(cdf, q, side="left")) for q in qs}
