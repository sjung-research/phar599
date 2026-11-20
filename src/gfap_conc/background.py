"""Analysis 3 - background fluorescence.

Without no-primary controls or annotated negative regions, background can only be
estimated by a **proxy**: tissue pixels that are below the shared GFAP foreground
threshold *and* lie at least ``exclusion_radius_px`` (in XY) away from any foreground
pixel. Every output row carries ``background_method`` so proxy values are never
mistaken for a measured autofluorescence level. Possible causes (perfusion,
non-specific binding, detector offset) cannot be distinguished by this measure.

Manual rectangular ROIs are given in **original pixel coordinates** as half-open
``[y0, y1) x [x0, x1)`` and are applied to every slice of the stack.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from scipy import ndimage as ndi
from skimage.morphology import disk

from .stats import robust_sd, select


def background_proxy_mask(
    foreground: np.ndarray,
    tissue_mask_2d: np.ndarray | None,
    exclusion_radius_px: int,
) -> np.ndarray:
    """Tissue & not within ``exclusion_radius_px`` (XY only) of foreground.

    Works for a 2-D image or slice-wise for a (Z, Y, X) foreground mask; the
    dilation never crosses slices, so z-sampling does not change the exclusion zone.
    """
    fg = np.asarray(foreground, dtype=bool)
    r = int(exclusion_radius_px)
    if fg.ndim == 2:
        near = ndi.binary_dilation(fg, structure=disk(r)) if r > 0 else fg
    elif fg.ndim == 3:
        near = ndi.binary_dilation(fg, structure=disk(r)[np.newaxis]) if r > 0 else fg
    else:
        raise ValueError(f"foreground mask must be 2-D or 3-D, got {fg.shape}")
    keep = ~near
    if tissue_mask_2d is not None:
        tm = np.asarray(tissue_mask_2d, dtype=bool)
        if tm.shape != fg.shape[-2:]:
            raise ValueError(f"tissue mask {tm.shape} does not match image YX {fg.shape[-2:]}")
        keep &= np.broadcast_to(tm, fg.shape)
    return keep


def roi_mask(shape_yx: tuple[int, int], rois: Sequence[Sequence[int]]) -> np.ndarray:
    """Union of half-open rectangles ``(y0, y1, x0, x1)`` in original pixel coordinates."""
    mask = np.zeros(shape_yx, dtype=bool)
    ny, nx = shape_yx
    for roi in rois:
        if len(roi) != 4:
            raise ValueError(f"ROI must be (y0, y1, x0, x1), got {roi}")
        y0, y1, x0, x1 = (int(v) for v in roi)
        if not (0 <= y0 < y1 <= ny and 0 <= x0 < x1 <= nx):
            raise ValueError(f"ROI {roi} lies outside the image ({ny} x {nx}); no implicit clipping")
        mask[y0:y1, x0:x1] = True
    return mask


def background_stats(values: np.ndarray, mask: np.ndarray, reference_count: int | None = None) -> dict:
    """Robust background summary of the masked pixels (raw units)."""
    v = select(values, mask)
    n = int(v.size)
    if n == 0:
        return {
            "bg_n_pixels": 0,
            "bg_fraction_of_reference": 0.0 if reference_count else np.nan,
            **dict.fromkeys(
                ["bg_median", "bg_mean", "bg_robust_sd", "bg_p05", "bg_p25", "bg_frac_zero"], np.nan
            ),
        }
    p05, p25 = np.percentile(v, [5.0, 25.0])
    return {
        "bg_n_pixels": n,
        "bg_fraction_of_reference": n / reference_count if reference_count else np.nan,
        "bg_median": float(np.median(v)),
        "bg_mean": float(np.mean(v, dtype=np.float64)),
        "bg_robust_sd": robust_sd(v),
        "bg_p05": float(p05),
        "bg_p25": float(p25),
        "bg_frac_zero": float(np.count_nonzero(v == 0)) / n,
    }


def signal_metrics(values: np.ndarray, foreground: np.ndarray, bg: dict) -> dict:
    """Foreground summary and (proxy-based) signal-to-background / contrast-to-noise."""
    fv = select(values, foreground)
    fg_median = float(np.median(fv)) if fv.size else np.nan
    fg_mean = float(np.mean(fv, dtype=np.float64)) if fv.size else np.nan
    bg_med, bg_sd = bg["bg_median"], bg["bg_robust_sd"]
    return {
        "fg_n_pixels": int(fv.size),
        "fg_median": fg_median,
        "fg_mean": fg_mean,
        "fg_median_minus_bg_median": fg_median - bg_med,
        "signal_to_background": fg_median / bg_med if bg_med and bg_med > 0 else np.nan,
        "contrast_to_noise": (fg_median - bg_med) / bg_sd if bg_sd and bg_sd > 0 else np.nan,
    }
