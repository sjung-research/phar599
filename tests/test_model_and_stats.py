"""Stack contract, immutability, and basic statistics."""

import numpy as np
import pytest
from conftest import make_stack

from gfap_conc.model import DeconvolvedIntensityError, ImageKind, StackValidationError, parse_dilution, require_raw
from gfap_conc.stats import robust_sd, summarize


def test_parse_dilution():
    assert parse_dilution("1:250") == 250
    assert parse_dilution(" 1:2000 ") == 2000
    with pytest.raises(ValueError):
        parse_dilution("250")


@pytest.mark.parametrize("shape", [(8, 8), (1, 2, 8, 8), (2, 1, 3, 8, 8)])
def test_non_3d_arrays_rejected(shape):
    """T/C axes must be resolved explicitly before a stack exists - no silent channel mixing."""
    with pytest.raises(StackValidationError):
        make_stack(np.zeros(shape, np.uint16))


def test_raw_must_be_integer_and_within_bit_depth():
    with pytest.raises(StackValidationError):
        make_stack(np.zeros((2, 4, 4), np.float32))
    with pytest.raises(StackValidationError, match="exceeds"):
        make_stack(np.full((2, 4, 4), 5000, np.uint16), bit_depth=12)
    # 16-bit container holding 12-bit data is fine
    s = make_stack(np.full((2, 4, 4), 4095, np.uint16), bit_depth=12)
    assert s.saturation_value == 4095


def test_stack_data_is_read_only_and_source_untouched():
    src = np.arange(2 * 4 * 4, dtype=np.uint16).reshape(2, 4, 4)
    s = make_stack(src)
    with pytest.raises(ValueError):
        s.data[0, 0, 0] = 1
    src_copy = src.copy()
    _ = summarize(s.data, s.saturation_value)
    np.testing.assert_array_equal(src, src_copy)
    assert src.flags.writeable  # caller's array is not modified


def test_require_raw_blocks_deconvolved():
    s = make_stack(np.zeros((2, 4, 4), np.uint16), kind=ImageKind.DECONVOLVED, bit_depth=16)
    with pytest.raises(DeconvolvedIntensityError):
        require_raw(s)


def test_summarize_known_values():
    v = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 4095], dtype=np.uint16)
    s = summarize(v, saturation_value=4095)
    assert s["n_pixels"] == 10
    assert s["mean"] == pytest.approx(v.astype(float).mean())
    assert s["median"] == 4.5
    assert s["variance"] == pytest.approx(v.astype(float).var())
    assert s["frac_zero"] == 0.1
    assert s["frac_saturated"] == 0.1
    assert s["p99_8"] == pytest.approx(np.percentile(v, 99.8))


def test_summarize_with_mask_and_shape_mismatch():
    img = np.arange(16, dtype=np.uint16).reshape(4, 4)
    m = np.zeros((4, 4), bool)
    m[0] = True
    assert summarize(img, 4095, m)["mean"] == pytest.approx(1.5)
    with pytest.raises(ValueError):
        summarize(img, 4095, np.ones((2, 2), bool))
    assert summarize(img, 4095, np.zeros((4, 4), bool))["n_pixels"] == 0


def test_single_outlier_moves_max_not_p998(rng):
    """Rose's reason for the 99.8th percentile: one crazy pixel must not dominate."""
    img = rng.integers(100, 200, size=(100, 100)).astype(np.uint16)
    before = summarize(img, 4095)
    img[50, 50] = 4095
    after = summarize(img, 4095)
    assert after["max"] - before["max"] > 3800
    assert abs(after["p99_8"] - before["p99_8"]) < 2


def test_no_overflow_on_large_uint16():
    img = np.full((50, 1000, 1000), 65535, dtype=np.uint16)
    s = summarize(img, 65535)
    assert s["mean"] == 65535.0
    assert s["frac_saturated"] == 1.0


def test_robust_sd_normal(rng):
    v = rng.normal(100, 10, 200_000)
    assert robust_sd(v) == pytest.approx(10, rel=0.02)
