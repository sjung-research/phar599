"""Tissue and GFAP-foreground masks.

Design rules:

* Every mask rule is fixed **before** looking at concentration-specific results and
  is applied identically to all concentrations.
* The default tissue mask is ``"none"`` (the whole field counts as tissue). At high
  magnification in spinal cord the field is often entirely inside the section,
  and an automatic threshold would then split tissue into "bright" and "dark"
  instead of tissue vs empty space. Switch to ``"absolute"``/``"otsu"`` only after
  inspecting the real images.
* The GFAP foreground threshold is a single **shared absolute** value (counts) for
  all concentrations. ``pooled_otsu`` derives that one value from an equal-sized
  pixel sample of every raw stack, so no concentration dominates it.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from scipy import ndimage as ndi
from skimage.filters import threshold_otsu
from skimage.morphology import remove_small_objects

TISSUE_METHODS = ("none", "absolute", "otsu")


def tissue_mask_2d(
    reference: np.ndarray,
    method: str = "none",
    threshold: float | None = None,
    sigma_px: float = 2.0,
    min_object_px: int = 500,
) -> np.ndarray:
    """2-D boolean tissue mask from a reference image (e.g. GFAP or DAPI MIP).

    The reference may be the DAPI channel: it is used only to locate tissue,
    never as a GFAP intensity measurement.
    """
    reference = np.asarray(reference)
    if reference.ndim != 2:
        raise ValueError(f"tissue reference must be 2-D, got {reference.shape}")
    if method == "none":
        return np.ones(reference.shape, dtype=bool)
    if method not in TISSUE_METHODS:
        raise ValueError(f"unknown tissue mask method {method!r}; choose from {TISSUE_METHODS}")
    smooth = ndi.gaussian_filter(reference.astype(np.float64), sigma_px)
    if method == "absolute":
        if threshold is None:
            raise ValueError("tissue mask method 'absolute' requires a threshold")
        level = float(threshold)
    else:
        level = float(threshold_otsu(smooth))
    mask = smooth >= level
    mask = remove_small_objects(mask, max_size=max(int(min_object_px) - 1, 0))
    return ndi.binary_fill_holes(mask)


def equal_subsample(arrays: Sequence[np.ndarray], n_per_array: int = 200_000, seed: int = 0) -> np.ndarray:
    """Deterministic, equal-sized random sample of pixel values from each array."""
    rng = np.random.default_rng(seed)
    parts = []
    for a in arrays:
        flat = np.asarray(a).ravel()
        n = min(n_per_array, flat.size)
        parts.append(flat[rng.choice(flat.size, size=n, replace=False)])
    return np.concatenate(parts)


def pooled_otsu_threshold(arrays: Sequence[np.ndarray], n_per_array: int = 200_000, seed: int = 0) -> float:
    """One Otsu threshold from pooled, equally weighted pixels of all raw stacks."""
    return float(threshold_otsu(equal_subsample(arrays, n_per_array, seed)))


def foreground_mask(data: np.ndarray, threshold: float) -> np.ndarray:
    """Pixels at or above a shared absolute threshold (raw counts)."""
    return np.asarray(data) >= threshold
