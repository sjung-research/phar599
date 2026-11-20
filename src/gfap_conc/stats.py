"""Descriptive pixel statistics shared by all analyses.

Percentiles use numpy's default ``method="linear"`` (type-7 interpolation); for
integer data the 99.8th percentile can therefore be non-integer. All accumulation
is done in float64 so that large uint16 stacks cannot overflow.
"""

from __future__ import annotations

import numpy as np

PERCENTILE_METHOD = "linear"

STAT_KEYS = (
    "n_pixels",
    "mean",
    "median",
    "std",
    "variance",
    "p95",
    "p99_8",
    "min",
    "max",
    "frac_zero",
    "frac_saturated",
)


def select(values: np.ndarray, mask: np.ndarray | None = None) -> np.ndarray:
    """Flattened values, optionally restricted to a boolean mask of identical shape."""
    values = np.asarray(values)
    if mask is None:
        return values.ravel()
    mask = np.asarray(mask, dtype=bool)
    if mask.shape != values.shape:
        raise ValueError(f"mask shape {mask.shape} does not match image shape {values.shape}")
    return values[mask]


def percentile(values: np.ndarray, q: float | list[float], mask: np.ndarray | None = None):
    v = select(values, mask)
    if v.size == 0:
        return np.nan if np.isscalar(q) else np.full(len(q), np.nan)
    return np.percentile(v, q, method=PERCENTILE_METHOD)


def summarize(
    values: np.ndarray,
    saturation_value: int | None,
    mask: np.ndarray | None = None,
) -> dict[str, float]:
    """Mean, median, spread, percentiles and clipping fractions of raw pixel values.

    ``variance`` is the population (ddof=0) variance of the selected pixels, i.e. a
    *within-image spatial* variance. It is not a variance across depth.
    ``frac_saturated`` counts pixels at or above ``saturation_value`` (2**bits - 1).
    """
    v = select(values, mask)
    n = int(v.size)
    if n == 0:
        out = dict.fromkeys(STAT_KEYS, np.nan)
        out["n_pixels"] = 0
        return out
    p95, p998 = np.percentile(v, [95.0, 99.8], method=PERCENTILE_METHOD)
    mean = float(np.mean(v, dtype=np.float64))
    var = float(np.var(v, dtype=np.float64))
    return {
        "n_pixels": n,
        "mean": mean,
        "median": float(np.median(v)),
        "std": float(np.sqrt(var)),
        "variance": var,
        "p95": float(p95),
        "p99_8": float(p998),
        "min": float(v.min()),
        "max": float(v.max()),
        "frac_zero": float(np.count_nonzero(v == 0)) / n,
        "frac_saturated": (
            float(np.count_nonzero(v >= saturation_value)) / n if saturation_value is not None else np.nan
        ),
    }


def robust_sd(values: np.ndarray) -> float:
    """Median absolute deviation scaled to the SD of a normal distribution (x1.4826)."""
    v = np.asarray(values).ravel()
    if v.size == 0:
        return np.nan
    med = np.median(v)
    return float(1.4826 * np.median(np.abs(v.astype(np.float64) - med)))
