"""Analysis 2: per-slice statistics, depth windows, middle-vs-outer retention."""

import numpy as np
import pytest
from conftest import make_stack

from gfap_conc.depth import depth_windows, detect_signal_z_extent, middle_vs_outer, slice_profile


def known_profile_stack(means, shape=(16, 16)):
    return make_stack(np.stack([np.full(shape, m, np.uint16) for m in means]), z_step_um=0.5)


def test_slice_profile_known_means_and_depth():
    s = known_profile_stack([100, 80, 50, 80, 120])
    p = slice_profile(s)
    assert list(p["mean"]) == [100, 80, 50, 80, 120]
    assert list(p["z_um"]) == [0, 0.5, 1.0, 1.5, 2.0]
    assert p["depth_fraction_stack"].iloc[-1] == 1.0
    assert (p["variance"] == 0).all()  # within-slice spatial variance of a constant slice
    assert (p["concentration"] == s.concentration).all()


def test_within_slice_variance_is_not_depth_variance():
    data = np.zeros((3, 2, 2), np.uint16)
    data[:, 0, :] = 10  # each slice: half 10, half 0 -> spatial variance 25
    s = make_stack(data)
    p = slice_profile(s)
    assert np.allclose(p["variance"], 25.0)
    r = middle_vs_outer(p, "mean", (0, 2), window_fraction=1 / 3)
    assert r["depth_cv"] == 0.0  # identical slices: no variation across depth


def test_foreground_fraction_uses_shared_threshold():
    data = np.zeros((2, 4, 4), np.uint16)
    data[0, :2] = 500
    p = slice_profile(make_stack(data), foreground_threshold=400)
    assert list(p["foreground_fraction"]) == [0.5, 0.0]


def test_depth_windows_non_overlapping():
    w = depth_windows(10, 0.2)
    assert list(w["first"]) == [0, 1] and list(w["last"]) == [8, 9] and list(w["middle"]) == [4, 5]
    with pytest.raises(ValueError):
        depth_windows(5, 0.3)  # k=2, needs >= 6 slices
    with pytest.raises(ValueError):
        depth_windows(10, 0.5)


def test_middle_vs_outer_asymmetric_profile():
    means = [200, 200, 100, 100, 50, 50, 100, 100, 150, 150]
    s = known_profile_stack(means)
    r = middle_vs_outer(slice_profile(s), "mean", (0, 9), 0.2)
    assert r["first_window_mean"] == 200
    assert r["middle_window_mean"] == 50
    assert r["last_window_mean"] == 150
    assert r["middle_over_first"] == pytest.approx(0.25)
    assert r["middle_over_last"] == pytest.approx(1 / 3)
    assert r["middle_over_outer_mean"] == pytest.approx(50 / 175)
    assert r["z_index_of_min"] == 4
    assert r["middle_window_slices"] == "4-5"


def test_middle_vs_outer_respects_z_range():
    means = [0, 0, 100, 100, 50, 50, 100, 100, 0, 0]
    r = middle_vs_outer(slice_profile(known_profile_stack(means)), "mean", (2, 7), 1 / 3)
    assert r["n_slices_analysed"] == 6
    assert r["first_window_mean"] == 100 and r["middle_window_mean"] == 50


def test_detect_signal_z_extent():
    v = np.array([1, 2, 50, 100, 90, 60, 3, 1], float)
    assert detect_signal_z_extent(v, 0.5) == (3, 5)
    assert detect_signal_z_extent(np.ones(4)) == (0, 3)
