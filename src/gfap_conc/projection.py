"""Analysis 1 - maximum-intensity projections (Rose: "average bit value on max-projected
images", with a 99.8th-percentile "normalization max").

Two quantities are deliberately kept in separate columns:

* ``mean`` (and the other stats) - computed on the *unscaled raw* MIP, in source units.
* ``display_vmax_within_image`` - the per-image 99.8th percentile, used only as the
  upper end of a display range. ``normalized_mean_secondary = mean / display_vmax``
  is reported as a secondary contrast/shape diagnostic; it removes absolute
  brightness differences by construction and must not be used to rank concentrations.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from .model import GFAPStack, require_raw
from .stats import percentile, summarize


def max_projection(data: np.ndarray) -> np.ndarray:
    """``MIP[y, x] = max_z I[z, y, x]`` for a (Z, Y, X) array; dtype is preserved."""
    data = np.asarray(data)
    if data.ndim != 3:
        raise ValueError(f"expected (Z, Y, X) array, got shape {data.shape}")
    return data.max(axis=0)


def window_projection(data: np.ndarray, z_start: int, z_stop: int) -> np.ndarray:
    """MIP over the half-open slice window ``[z_start, z_stop)``."""
    n_z = np.asarray(data).shape[0]
    if not 0 <= z_start < z_stop <= n_z:
        raise ValueError(f"invalid z window [{z_start}, {z_stop}) for {n_z} slices")
    return max_projection(np.asarray(data)[z_start:z_stop])


def matched_depth_windows(
    stacks: Sequence[GFAPStack],
    z_ranges: dict[str, tuple[int, int]],
) -> dict[str, dict[str, tuple[int, int]]]:
    """Projection windows covering the *same physical thickness* in every stack.

    The thickness is that of the thinnest analysed z-range. Two placements are
    returned per stack: ``"first"`` (starting at the first analysed slice, i.e. the
    first acquired side) and ``"center"`` (centred in the analysed range).
    Windows are half-open ``(start, stop)`` slice indices. If any stack lacks a
    z-step, slice counts are matched instead of micrometres.
    """
    use_um = all(s.z_step_um is not None for s in stacks)
    extents = {}
    for s in stacks:
        z0, z1 = z_ranges[s.series_name]
        n = z1 - z0 + 1
        extents[s.series_name] = (n - 1) * s.z_step_um if use_um else float(n - 1)
    target = min(extents.values())

    out: dict[str, dict[str, tuple[int, int]]] = {}
    for s in stacks:
        z0, z1 = z_ranges[s.series_name]
        n_avail = z1 - z0 + 1
        step = s.z_step_um if use_um else 1.0
        n_win = int(np.floor(target / step + 1e-9)) + 1
        n_win = max(1, min(n_win, n_avail))
        c0 = z0 + (n_avail - n_win) // 2
        out[s.series_name] = {"first": (z0, z0 + n_win), "center": (c0, c0 + n_win)}
    return out


def projection_rows(
    stack: GFAPStack,
    mip: np.ndarray,
    projection: str,
    z_window: tuple[int, int],
    tissue_mask: np.ndarray | None,
    display_percentile: float = 99.8,
) -> list[dict]:
    """Summary rows (region = all pixels, and tissue if a mask is supplied)."""
    require_raw(stack, "MIP intensity summaries")
    regions: list[tuple[str, np.ndarray | None]] = [("all_pixels", None)]
    if tissue_mask is not None and not tissue_mask.all():
        regions.append(("tissue_mask", tissue_mask))
    rows = []
    for region, mask in regions:
        s = summarize(mip, stack.saturation_value, mask)
        vmax = float(percentile(mip, display_percentile, mask))
        rows.append(
            {
                **stack.identity(),
                "projection": projection,
                "z_start": int(z_window[0]),
                "z_stop_exclusive": int(z_window[1]),
                "n_z_projected": int(z_window[1] - z_window[0]),
                "z_thickness_um": (
                    (z_window[1] - z_window[0] - 1) * stack.z_step_um if stack.z_step_um else np.nan
                ),
                "region": region,
                **s,
                "display_percentile": display_percentile,
                "display_vmax_within_image": vmax,
                "normalized_mean_secondary": s["mean"] / vmax if vmax > 0 else np.nan,
            }
        )
    return rows


def shared_display_vmax(mips: Sequence[np.ndarray], display_percentile: float = 99.8) -> float:
    """One display ceiling for all panels: the largest per-image percentile.

    Using the maximum keeps the brightest image unsaturated on screen and shows
    dimmer images as dimmer - the honest cross-concentration view.
    """
    return float(max(percentile(m, display_percentile) for m in mips))
