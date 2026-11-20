"""Analysis 1: MIP axis, raw vs normalised separation, matched depth windows."""

import numpy as np
import pytest
from conftest import make_stack

from gfap_conc.model import DeconvolvedIntensityError, ImageKind
from gfap_conc.projection import (
    matched_depth_windows,
    max_projection,
    projection_rows,
    shared_display_vmax,
    window_projection,
)


def test_mip_is_max_over_z_axis_only(rng):
    data = rng.integers(0, 4000, size=(5, 6, 7)).astype(np.uint16)
    mip = max_projection(data)
    assert mip.shape == (6, 7)
    assert mip.dtype == data.dtype
    np.testing.assert_array_equal(mip, data.max(axis=0))
    # a non-cubic shape makes the wrong axis impossible to match by accident
    assert not np.array_equal(mip.shape, data.max(axis=1).shape)


def test_window_projection_bounds():
    data = np.arange(4 * 2 * 2, dtype=np.uint16).reshape(4, 2, 2)
    np.testing.assert_array_equal(window_projection(data, 1, 3), data[2])
    with pytest.raises(ValueError):
        window_projection(data, 2, 2)
    with pytest.raises(ValueError):
        window_projection(data, 0, 5)


def test_per_image_normalisation_hides_twofold_difference(rng):
    """A stack that is exactly 2x brighter must stay 2x brighter in the raw columns,
    while the secondary normalised mean is (as designed) identical."""
    base = rng.integers(50, 1000, size=(6, 40, 40)).astype(np.uint16)
    dim = make_stack(base, "1:2000")
    bright = make_stack(base * 2, "1:250")
    r_dim = projection_rows(dim, max_projection(dim.data), "mip", (0, 6), None)[0]
    r_bright = projection_rows(bright, max_projection(bright.data), "mip", (0, 6), None)[0]
    assert r_bright["mean"] == pytest.approx(2 * r_dim["mean"])
    assert r_bright["display_vmax_within_image"] == pytest.approx(2 * r_dim["display_vmax_within_image"])
    assert r_bright["normalized_mean_secondary"] == pytest.approx(r_dim["normalized_mean_secondary"])
    # raw and display columns are distinct quantities
    assert r_dim["mean"] != r_dim["display_vmax_within_image"]


def test_projection_rows_carry_identity_and_reject_deconvolved():
    s = make_stack(np.ones((2, 4, 4), np.uint16), "1:1000", name="Series 7")
    row = projection_rows(s, max_projection(s.data), "mip_full_stack", (0, 2), None)[0]
    assert (row["concentration"], row["series_name"], row["channel_index"], row["image_kind"]) == (
        "1:1000", "Series 7", 0, "raw")
    d = make_stack(np.ones((2, 4, 4), np.uint16), "1:1000", kind=ImageKind.DECONVOLVED, bit_depth=16)
    with pytest.raises(DeconvolvedIntensityError):
        projection_rows(d, max_projection(d.data), "mip", (0, 2), None)


def test_tissue_region_row_added_only_for_partial_mask():
    s = make_stack(np.arange(32, dtype=np.uint16).reshape(2, 4, 4))
    mip = max_projection(s.data)
    assert len(projection_rows(s, mip, "mip", (0, 2), np.ones((4, 4), bool))) == 1
    m = np.zeros((4, 4), bool)
    m[:2] = True
    rows = projection_rows(s, mip, "mip", (0, 2), m)
    assert [r["region"] for r in rows] == ["all_pixels", "tissue_mask"]
    assert rows[1]["n_pixels"] == 8


def test_shared_display_vmax_is_largest_percentile():
    a = np.full((4, 4), 100, np.uint16)
    b = np.full((4, 4), 300, np.uint16)
    assert shared_display_vmax([a, b]) == 300


def test_matched_depth_windows_physical_thickness():
    thin = make_stack(np.zeros((11, 4, 4), np.uint16), "1:250", z_step_um=1.0)  # 10 µm
    thick = make_stack(np.zeros((41, 4, 4), np.uint16), "1:500", z_step_um=0.5)  # 20 µm
    w = matched_depth_windows([thin, thick], {thin.series_name: (0, 10), thick.series_name: (0, 40)})
    assert w[thin.series_name]["first"] == (0, 11)
    assert w[thick.series_name]["first"] == (0, 21)  # 21 slices x 0.5 µm = 10 µm span
    c0, c1 = w[thick.series_name]["center"]
    assert c1 - c0 == 21 and c0 == 10
