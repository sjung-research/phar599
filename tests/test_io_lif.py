"""Reader adapter on synthetic LIF files written to the documented layout."""

import numpy as np
import pytest
from lif_writer import LifTestImage, write_lif

from gfap_conc import io_lif
from gfap_conc.io_lif import channel_zyx, guess_concentration, guess_kind
from gfap_conc.model import ImageKind


@pytest.fixture
def lif_path(tmp_path, rng):
    images = []
    for conc in ("250", "500"):
        d = rng.integers(0, 4096, size=(6, 2, 16, 20)).astype(np.uint16)  # Z, C, Y, X
        images.append(LifTestImage(f"GFAP 1-{conc}", d, hardware_attrs={"Zoom": "2", "Pinhole": "1.0e-4"},
                                   detectors=[{"Name": "HyD 1", "IsActive": "1", "Gain": "100"},
                                              {"Name": "PMT 2", "IsActive": "0", "Gain": "700"}]))
        images.append(LifTestImage(f"GFAP 1-{conc}_lng", d, resolution_bits=16))
    images[2].hardware_attrs = {"Zoom": "3", "Pinhole": "1.0e-4"}  # 1:500 raw differs in zoom
    p = write_lif(tmp_path / "curve.lif", images)
    return p, images


def test_list_series_axes_and_physical_sizes(lif_path):
    p, images = lif_path
    metas = io_lif.list_series(p)
    assert [m.path for m in metas] == [im.name for im in images]
    m = metas[0]
    assert m.sizes == {"Z": 6, "C": 2, "Y": 16, "X": 20}
    assert m.z_step_um == pytest.approx(0.5)
    assert m.pixel_size_um == pytest.approx((0.2, 0.2))
    assert [c.resolution_bits for c in m.channels] == [12, 12]
    assert [c.lut_name for c in m.channels] == ["Red", "Blue"]
    assert m.hardware["ATLConfocalSettingDefinition[0].Zoom"] == "2"
    assert "Detector[HyD 1].Gain" in m.hardware
    assert not any(k.startswith("Detector[PMT 2]") for k in m.hardware)  # inactive skipped


def test_load_stack_selects_channel_without_mixing(lif_path):
    p, images = lif_path
    for ch in (0, 1):
        st = io_lif.load_stack(p, "GFAP 1-250", ch, "1:250", "raw")
        np.testing.assert_array_equal(st.data, images[0].data[:, ch])
        assert st.data.shape == (6, 16, 20)
        assert st.bit_depth == 12 and st.saturation_value == 4095
        assert st.source.endswith("::GFAP 1-250")


def test_inventory_and_comparability(lif_path):
    p, images = lif_path
    inv, settings = io_lif.inventory(p)
    assert len(inv) == 8  # 4 series x 2 channels
    row = inv[(inv.series_path == "GFAP 1-250") & (inv.channel_index == 0)].iloc[0]
    assert row.raw_max == images[0].data[:, 0].max()
    assert row.image_kind_guess == "raw" and row.concentration_guess == "1:250"
    assert inv[inv.series_path == "GFAP 1-250_lng"].image_kind_guess.iloc[0] == "deconvolved"
    comp = io_lif.comparability_table(settings, ["GFAP 1-250", "GFAP 1-500"])
    flagged = set(comp.loc[comp.differs, "key"])
    assert "ATLConfocalSettingDefinition[0].Zoom" in flagged
    assert "ATLConfocalSettingDefinition[0].Pinhole" not in flagged


def test_sha256_stable(lif_path):
    p, _ = lif_path
    assert io_lif.sha256_file(p) == io_lif.sha256_file(p)
    assert len(io_lif.sha256_file(p)) == 64


def test_channel_zyx_rejects_unsupported_axes():
    with pytest.raises(ValueError, match="unsupported"):
        channel_zyx(np.zeros((2, 3, 4, 5)), ("T", "Z", "Y", "X"), 0)
    with pytest.raises(ValueError, match="no Z"):
        channel_zyx(np.zeros((2, 4, 5)), ("C", "Y", "X"), 0)
    with pytest.raises(ValueError, match="out of range"):
        channel_zyx(np.zeros((3, 2, 4, 5)), ("Z", "C", "Y", "X"), 2)
    a = np.arange(2 * 3 * 4 * 5).reshape(2, 3, 4, 5)  # C, Z, Y, X order is handled too
    np.testing.assert_array_equal(channel_zyx(a, ("C", "Z", "Y", "X"), 1), a[1])


@pytest.mark.parametrize("name,expected", [
    ("Ms GFAP 1:250", "1:250"), ("GFAP_1-2000_lng", "1:2000"), ("series 500", "1:500"),
    ("Position 3", None), ("1:250 vs 1:500", None), ("GFAP 1_1000 z", "1:1000"),
])
def test_guess_concentration(name, expected):
    assert guess_concentration(name) == expected


def test_guess_kind_uses_exact_token():
    assert guess_kind("A_lng") is ImageKind.DECONVOLVED
    assert guess_kind("A_LNG") is ImageKind.RAW  # not normalised; inventory flags it for review
