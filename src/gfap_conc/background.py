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
                ["bg_median", "bg_mean", "bg_robust_sd", "bg_p05", "bg_p25", "bg_p95", "bg_tail_asymmetry",
                 "bg_frac_zero"], np.nan
            ),
        }
    p05, p25, p50, p95 = np.percentile(v, [5.0, 25.0, 50.0, 95.0])
    return {
        "bg_n_pixels": n,
        "bg_fraction_of_reference": n / reference_count if reference_count else np.nan,
        "bg_median": float(p50),
        "bg_mean": float(np.mean(v, dtype=np.float64)),
        "bg_robust_sd": robust_sd(v),
        "bg_p05": float(p05),
        "bg_p25": float(p25),
        "bg_p95": float(p95),
        "bg_tail_asymmetry": tail_asymmetry(p05, p50, p95),
        "bg_frac_zero": float(np.count_nonzero(v == 0)) / n,
    }


def tail_asymmetry(p05: float, p50: float, p95: float) -> float:
    """``(p95 - median) / (median - p05)``: ~1 for symmetric noise, >1 for a right tail.

    A pure background (detector noise, uniform autofluorescence) is roughly
    symmetric; faint GFAP processes that fall below the shared threshold add a
    right tail to the proxy distribution.
    """
    lower = p50 - p05
    return float((p95 - p50) / lower) if lower > 0 else np.nan


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


PROXY_NOTE = ("proxy: tissue voxels below the shared foreground threshold and away from foreground; "
              "NOT a measurement of autofluorescence (cause cannot be determined from these images)")
ROI_NOTE = "manual ROI in original pixel coordinates; content assumed GFAP-free by the user, not verified"


def threshold_sensitivity_rows(
    data: np.ndarray,
    tissue_mask_2d: np.ndarray | None,
    thresholds: Sequence[tuple[float, float]],
    exclusion_radius_px: int,
) -> list[dict]:
    """Re-estimate foreground and background proxy at several shared thresholds.

    ``thresholds`` are ``(factor, absolute_threshold)`` pairs. Diagnostic only: the
    analysis itself always uses factor 1.0. A background median that stays flat across
    reasonable thresholds supports the proxy; one that climbs with the threshold means
    sub-threshold GFAP signal is being counted as background.
    """
    data = np.asarray(data)
    tm = None if tissue_mask_2d is None else np.asarray(tissue_mask_2d, bool)
    tissue3d = np.ones(data.shape, bool) if tm is None else np.broadcast_to(tm, data.shape)
    n_tissue = int(np.count_nonzero(tissue3d))
    rows = []
    for factor, thr in thresholds:
        fg = (data >= thr) & tissue3d
        bg_mask = background_proxy_mask(fg, tm, exclusion_radius_px)
        bg = background_stats(data, bg_mask, n_tissue)
        rows.append({
            "threshold_factor": float(factor),
            "foreground_threshold": float(thr),
            "fg_fraction_of_tissue": int(np.count_nonzero(fg)) / n_tissue if n_tissue else np.nan,
            **bg,
            **signal_metrics(data, fg, bg),
        })
    return rows


def threshold_dependence(sens_rows: Sequence[dict]) -> float:
    """Relative change of the proxy median between the factors bracketing 1.0.

    ``(median(f_hi) - median(f_lo)) / median(1.0)`` with f_lo/f_hi the nearest factors
    below/above 1.0. NaN if either side is missing.
    """
    by_f = {r["threshold_factor"]: r["bg_median"] for r in sens_rows}
    lo = [f for f in by_f if f < 1.0]
    hi = [f for f in by_f if f > 1.0]
    ref = by_f.get(1.0, np.nan)
    if not lo or not hi or not ref or not np.isfinite(ref) or ref <= 0:
        return np.nan
    return float((by_f[min(hi)] - by_f[max(lo)]) / ref)


def background_qc_flags(
    bg: dict,
    dependence: float,
    max_fraction: float,
    max_tail_asymmetry: float,
    max_dependence: float,
) -> list[str]:
    """Machine-readable reasons to distrust a background-proxy estimate."""
    flags = []
    if bg.get("bg_n_pixels", 0) == 0:
        return ["EMPTY_REGION"]
    frac = bg.get("bg_fraction_of_reference", np.nan)
    if np.isfinite(frac) and frac > max_fraction:
        flags.append("LARGE_REGION")
    if np.isfinite(frac) and frac < 0.05:
        flags.append("SMALL_REGION")
    tail = bg.get("bg_tail_asymmetry", np.nan)
    if np.isfinite(tail) and tail > max_tail_asymmetry:
        flags.append("RIGHT_TAIL")
    if np.isfinite(dependence) and abs(dependence) > max_dependence:
        flags.append("THRESHOLD_DEPENDENT")
    return flags


FLAG_MEANING = {
    "EMPTY_REGION": "no voxels qualify as background; background and SBR undefined",
    "SMALL_REGION": "fewer than 5% of tissue voxels qualify; estimate rests on few, possibly atypical voxels",
    "LARGE_REGION": "most tissue voxels are classed as background; the shared threshold may lie above much "
                    "of this stain's real GFAP signal",
    "RIGHT_TAIL": "proxy intensities have a long upper tail, as expected if faint GFAP processes are included",
    "THRESHOLD_DEPENDENT": "proxy median changes with the threshold, i.e. it contains sub-threshold signal",
}
