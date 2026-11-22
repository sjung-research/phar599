"""Analysis 2 - intensity changes through the z-stack ("how much darker is it in the middle?").

Definitions (also in docs/methods.md):

* Per-slice statistics are computed on raw values of each z-slice (optionally within
  the 2-D tissue mask). ``variance`` is the within-slice *spatial* pixel variance.
* Depth variation is summarised separately from the per-slice series (``depth_cv``,
  ``min_over_max``); the two notions of "variance" are never mixed.
* z = 0 is the **first acquired slice**. Whether that is the coverslip side is not
  known from pixel data, so windows are named ``first`` / ``last``, not top/bottom.
* Middle-vs-outer windows: within the analysed range of n slices,
  ``k = max(1, round(n * window_fraction))``; first = slices [0, k), last =
  [n-k, n), middle = k slices centred at (n-k)//2. Requires n >= 3k so the
  windows cannot overlap.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .model import GFAPStack, require_raw
from .stats import summarize


def slice_profile(
    stack: GFAPStack,
    tissue_mask: np.ndarray | None = None,
    foreground_threshold: float | None = None,
    z_range: tuple[int, int] | None = None,
) -> pd.DataFrame:
    """One row per z-slice with raw statistics and depth coordinates."""
    require_raw(stack, "depth profiles")
    n_z = stack.n_z
    z0, z1 = z_range if z_range is not None else (0, n_z - 1)
    z_um = stack.z_um()
    rows = []
    for z in range(n_z):
        plane = stack.data[z]
        s = summarize(plane, stack.saturation_value, tissue_mask)
        in_range = z0 <= z <= z1
        row = {
            **stack.identity(),
            "z_index": z,
            "z_um": float(z_um[z]),
            "depth_fraction_stack": z / (n_z - 1) if n_z > 1 else 0.0,
            "in_analysis_range": in_range,
            "depth_fraction_analysis": ((z - z0) / (z1 - z0) if z1 > z0 else 0.0) if in_range else np.nan,
            "region": "all_pixels" if tissue_mask is None or tissue_mask.all() else "tissue_mask",
            **s,
        }
        if foreground_threshold is not None:
            v = plane if tissue_mask is None else plane[tissue_mask]
            row["foreground_threshold"] = float(foreground_threshold)
            row["foreground_fraction"] = float(np.count_nonzero(v >= foreground_threshold)) / max(v.size, 1)
        rows.append(row)
    return pd.DataFrame(rows)


def detect_signal_z_extent(
    values: np.ndarray, rel_level: float = 0.2, floor: float | None = None
) -> tuple[int, int]:
    """First and last slice whose value reaches ``floor + rel_level * (max - floor)``.

    Meant to trim slices acquired outside the tissue (e.g. in mounting medium). Only
    the ends can be trimmed: every slice between the first and last crossing is kept,
    so a mid-stack intensity dip is never removed. ``floor`` should be a background
    level (the pipeline passes the lowest per-slice median); if omitted, the profile
    minimum is used, which is only appropriate for profiles without a deep internal dip.
    Returned indices are inclusive. **Unvalidated on real data** - review the depth
    profile before enabling it.

    History: with ``floor = min(profile)`` and ``rel_level = 0.5``, the rule cut the real
    stacks in the middle (their two bright surfaces and mid-stack dip put the half level
    above the second peak); hence the background floor and the lower default level.
    """
    v = np.asarray(values, dtype=np.float64)
    if v.ndim != 1 or v.size == 0:
        raise ValueError("expected a non-empty 1-D profile")
    hi = np.nanmax(v)
    lo = np.nanmin(v) if floor is None else float(floor)
    if not hi > lo:
        return 0, v.size - 1
    idx = np.flatnonzero(v >= lo + rel_level * (hi - lo))
    return int(idx[0]), int(idx[-1])


def depth_windows(n: int, window_fraction: float = 0.2) -> dict[str, np.ndarray]:
    """Relative slice indices of the first, middle and last windows (see module doc)."""
    if not 0 < window_fraction < 1 / 3 + 1e-9:
        raise ValueError("window_fraction must be in (0, 1/3]")
    k = max(1, int(round(n * window_fraction)))
    if 3 * k > n:
        raise ValueError(f"{n} slices are too few for three non-overlapping windows of {k}")
    m0 = (n - k) // 2
    return {
        "first": np.arange(0, k),
        "middle": np.arange(m0, m0 + k),
        "last": np.arange(n - k, n),
    }


def middle_vs_outer(
    profile: pd.DataFrame,
    metric: str,
    z_range: tuple[int, int],
    window_fraction: float = 0.2,
) -> dict:
    """Signal retention in the middle of the stack relative to each outer window.

    Ratios use the mean of the per-slice ``metric`` within each window. The first
    and last windows are reported separately because profiles are rarely symmetric.
    """
    z0, z1 = z_range
    sub = profile[(profile["z_index"] >= z0) & (profile["z_index"] <= z1)].sort_values("z_index")
    vals = sub[metric].to_numpy(dtype=np.float64)
    zs = sub["z_index"].to_numpy()
    w = depth_windows(len(vals), window_fraction)
    first, middle, last = (float(np.mean(vals[w[k]])) for k in ("first", "middle", "last"))
    outer = (first + last) / 2

    def ratio(a: float, b: float) -> float:
        return a / b if b > 0 else np.nan

    def span(key: str) -> str:
        idx = zs[w[key]]
        return f"{idx[0]}-{idx[-1]}"

    mean_all = float(np.mean(vals))
    return {
        "metric": metric,
        "z_first_analysed": int(z0),
        "z_last_analysed": int(z1),
        "n_slices_analysed": int(len(vals)),
        "window_fraction": window_fraction,
        "window_n_slices": int(len(w["first"])),
        "first_window_slices": span("first"),
        "middle_window_slices": span("middle"),
        "last_window_slices": span("last"),
        "first_window_mean": first,
        "middle_window_mean": middle,
        "last_window_mean": last,
        "middle_over_first": ratio(middle, first),
        "middle_over_last": ratio(middle, last),
        "middle_over_outer_mean": ratio(middle, outer),
        "depth_cv": float(np.std(vals) / mean_all) if mean_all > 0 else np.nan,
        "min_over_max": ratio(float(np.min(vals)), float(np.max(vals))),
        "z_index_of_min": int(zs[int(np.argmin(vals))]),
        "z_index_of_max": int(zs[int(np.argmax(vals))]),
    }
