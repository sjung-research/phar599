"""Analysis 3 + masks: exclusion zones, ROI coordinates, shared thresholds."""

import numpy as np
import pytest

from gfap_conc.background import background_proxy_mask, background_stats, roi_mask, signal_metrics
from gfap_conc.masks import equal_subsample, foreground_mask, pooled_otsu_threshold, tissue_mask_2d


def test_background_proxy_excludes_dilated_foreground_xy_only():
    fg = np.zeros((3, 21, 21), bool)
    fg[1, 10, 10] = True
    bg = background_proxy_mask(fg, None, exclusion_radius_px=3)
    assert not bg[1, 10, 13] and bg[1, 10, 14]  # radius 3 in XY
    assert bg[0, 10, 10] and bg[2, 10, 10]  # no dilation across slices


def test_background_proxy_respects_tissue_mask():
    fg = np.zeros((10, 10), bool)
    tissue = np.zeros((10, 10), bool)
    tissue[:, :5] = True
    bg = background_proxy_mask(fg, tissue, 0)
    assert bg[:, :5].all() and not bg[:, 5:].any()
    with pytest.raises(ValueError):
        background_proxy_mask(fg, np.ones((5, 5), bool), 1)


def test_roi_mask_uses_original_pixel_coordinates():
    img = np.arange(100).reshape(10, 10)
    m = roi_mask(img.shape, [(2, 4, 5, 8)])
    assert m.sum() == 6
    np.testing.assert_array_equal(img[m], [25, 26, 27, 35, 36, 37])
    with pytest.raises(ValueError, match="outside"):
        roi_mask(img.shape, [(0, 11, 0, 2)])  # no silent clipping / resizing


def test_background_stats_and_snr():
    img = np.full((10, 10), 50, np.uint16)
    img[:, :2] = 500
    fg = img >= 400
    bg = background_stats(img, ~fg, reference_count=100)
    assert bg["bg_median"] == 50 and bg["bg_robust_sd"] == 0 and bg["bg_fraction_of_reference"] == 0.8
    sm = signal_metrics(img, fg, bg)
    assert sm["signal_to_background"] == 10
    assert np.isnan(sm["contrast_to_noise"])  # zero spread: undefined, not infinite
    empty = background_stats(img, np.zeros_like(fg))
    assert empty["bg_n_pixels"] == 0 and np.isnan(empty["bg_median"])


def test_tissue_mask_methods():
    img = np.zeros((60, 60))
    img[10:50, 10:50] = 100
    assert tissue_mask_2d(img, "none").all()
    m = tissue_mask_2d(img, "absolute", threshold=50, sigma_px=1, min_object_px=10)
    assert m[30, 30] and not m[2, 2]
    assert tissue_mask_2d(img, "otsu", sigma_px=1, min_object_px=10)[30, 30]
    with pytest.raises(ValueError):
        tissue_mask_2d(img, "absolute")


def test_pooled_threshold_is_shared_and_equally_weighted(rng):
    small = rng.normal(100, 5, 1000)
    large = rng.normal(1000, 5, 100_000)
    sample = equal_subsample([small, large], n_per_array=1000)
    assert sample.size == 2000  # the bigger array does not dominate
    t = pooled_otsu_threshold([small, large], n_per_array=1000)
    assert 100 < t < 1000
    assert t == pooled_otsu_threshold([small, large], n_per_array=1000)  # deterministic
    assert foreground_mask(np.array([t - 1, t]), t).tolist() == [False, True]
