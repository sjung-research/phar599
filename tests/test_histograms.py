"""Analysis 4: shared integer-aligned bins, exact totals, survival curves."""

import numpy as np
import pytest

from gfap_conc.histograms import (
    binned_counts,
    histogram_table,
    integer_bin_edges,
    level_counts,
    quantiles_from_levels,
    survival_table,
)


@pytest.mark.parametrize("max_value,n_bins", [(4095, 256), (65535, 1024), (255, 256), (4095, 300)])
def test_integer_bin_edges(max_value, n_bins):
    e = integer_bin_edges(max_value, n_bins)
    assert len(e) <= n_bins + 1 and e[0] == 0 and e[-1] >= max_value + 1
    assert e[-2] <= max_value  # no bin lies entirely above the representable range
    assert np.all(np.diff(e) == e[1] - e[0])  # equal integer width
    assert e.dtype.kind == "i"


def test_counts_match_numpy_histogram_and_totals(rng):
    v = rng.integers(0, 4096, size=(7, 33, 31)).astype(np.uint16)
    e = integer_bin_edges(4095, 256)
    levels = level_counts(v, 4095)
    counts = binned_counts(levels, e)
    np.testing.assert_array_equal(counts, np.histogram(v, bins=e)[0])
    t = histogram_table(levels, e, 4095)
    assert t["count"].sum() == v.size
    assert t["fraction"].sum() == pytest.approx(1.0)
    assert (t["density_per_level"] * t["n_levels_in_bin"]).sum() == pytest.approx(1.0)


def test_last_partial_bin_density_corrected():
    e = integer_bin_edges(4095, 300)  # width 14 -> last bin holds fewer valid levels
    t = histogram_table(np.ones(4096, np.int64), e, 4095)
    assert t["n_levels_in_bin"].iloc[-1] < t["n_levels_in_bin"].iloc[0]
    assert np.allclose(t["density_per_level"], 1 / 4096)  # uniform data -> flat density


def test_same_edges_for_different_conditions(rng):
    """Two concentrations with different ranges still share identical bins."""
    e = integer_bin_edges(4095, 128)
    a = histogram_table(level_counts(rng.integers(0, 300, 1000).astype(np.uint16), 4095), e, 4095)
    b = histogram_table(level_counts(rng.integers(0, 4000, 1000).astype(np.uint16), 4095), e, 4095)
    np.testing.assert_array_equal(a["bin_left"], b["bin_left"])


def test_level_counts_rejects_bad_input():
    with pytest.raises(TypeError):
        level_counts(np.array([1.5, 2.0]), 10)
    with pytest.raises(ValueError):
        level_counts(np.array([0, 11], np.uint16), 10)


def test_masked_histogram_total():
    img = np.arange(16, dtype=np.uint16).reshape(4, 4)
    m = img % 2 == 0
    assert level_counts(img, 4095, m).sum() == 8


def test_survival_and_quantiles():
    levels = level_counts(np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9], np.uint16), 15)
    s = survival_table(levels, np.array([0, 5, 9, 10]))
    assert s["fraction_at_or_above"].tolist() == [1.0, 0.5, 0.1, 0.0]
    q = quantiles_from_levels(levels, (0.5, 0.9))
    assert q == {"q0.5": 4, "q0.9": 8}
