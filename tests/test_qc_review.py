"""Background-proxy QC, threshold sensitivity, review hashing, origin separation, plot helpers."""

import json

import numpy as np
import pytest
from conftest import make_stack

from gfap_conc.background import (
    background_qc_flags,
    background_stats,
    tail_asymmetry,
    threshold_dependence,
    threshold_sensitivity_rows,
)
from gfap_conc.cli import main
from gfap_conc.config import AnalysisSettings, DepthSettings, ProjectConfig, SeriesSpec, TissueMaskSettings
from gfap_conc.pipeline import data_origin, resolve_z_ranges, run_analysis
from gfap_conc.plots import element_word, spread_positions
from gfap_conc.review import review_settings_sha256, review_status
from gfap_conc.synthetic import make_concentration_series


@pytest.fixture(scope="module")
def synthetic_result():
    return run_analysis(make_concentration_series(shape=(15, 64, 64)), AnalysisSettings())


# ----------------------------------------------------------- tail asymmetry / flags
def test_tail_asymmetry_symmetric_vs_contaminated(rng):
    noise = rng.poisson(60, 200_000)
    s = background_stats(noise, np.ones(noise.shape, bool))
    assert 0.7 < s["bg_tail_asymmetry"] < 1.4
    contaminated = np.concatenate([noise, rng.uniform(100, 300, 40_000)])
    c = background_stats(contaminated, np.ones(contaminated.shape, bool))
    assert c["bg_tail_asymmetry"] > 2.0
    assert np.isnan(tail_asymmetry(5, 5, 9))


def test_qc_flags_each_condition():
    base = {"bg_n_pixels": 100, "bg_fraction_of_reference": 0.5, "bg_tail_asymmetry": 1.0}
    assert background_qc_flags(base, 0.0, 0.85, 2.0, 0.1) == []
    assert background_qc_flags({**base, "bg_fraction_of_reference": 0.9}, 0.0, 0.85, 2.0, 0.1) == ["LARGE_REGION"]
    assert background_qc_flags({**base, "bg_fraction_of_reference": 0.01}, 0.0, 0.85, 2.0, 0.1) == ["SMALL_REGION"]
    assert background_qc_flags({**base, "bg_tail_asymmetry": 3.0}, 0.0, 0.85, 2.0, 0.1) == ["RIGHT_TAIL"]
    assert background_qc_flags(base, 0.25, 0.85, 2.0, 0.1) == ["THRESHOLD_DEPENDENT"]
    assert background_qc_flags({"bg_n_pixels": 0}, np.nan, 0.85, 2.0, 0.1) == ["EMPTY_REGION"]


# ------------------------------------------------------------- threshold sensitivity
def test_sensitivity_factor_one_matches_main_background_row(synthetic_result):
    t = synthetic_result.tables
    sens = t["threshold_sensitivity"]
    main_bg = t["background_summary"]
    main_bg = main_bg[main_bg.background_method == "proxy_stack_3d"].set_index("series_name")
    one = sens[sens.threshold_factor == 1.0].set_index("series_name")
    for col in ("bg_median", "bg_n_pixels", "bg_robust_sd", "fg_median", "signal_to_background"):
        np.testing.assert_allclose(one.loc[main_bg.index, col], main_bg[col])
    assert set(sens.groupby("series_name")["threshold_factor"].apply(tuple)) == {(0.5, 0.75, 1.0, 1.25, 1.5, 2.0)}
    # thresholds are shared absolute values: identical for every concentration at each factor
    assert sens.groupby("threshold_factor")["foreground_threshold"].nunique().eq(1).all()


def test_sensitivity_monotone_region_and_dependence():
    data = np.full((4, 40, 40), 50, np.uint16)
    data[:, 5:8, :] = 400   # bright process
    data[:, 20:22, :] = 150  # faint process: becomes 'background' once T > 150
    rows = threshold_sensitivity_rows(data, None, [(0.5, 100), (1.0, 200), (2.0, 400)], exclusion_radius_px=2)
    frac = [r["bg_fraction_of_reference"] for r in rows]
    assert frac == sorted(frac)  # higher threshold -> more voxels classed as background
    assert rows[0]["bg_median"] == 50
    dep = threshold_dependence(rows)
    assert np.isfinite(dep)
    assert threshold_dependence([{"threshold_factor": 1.0, "bg_median": 5.0}]) != threshold_dependence(rows)
    assert np.isnan(threshold_dependence([{"threshold_factor": 1.0, "bg_median": 5.0}]))


def test_synthetic_weak_stain_background_flagged_not_called_autofluorescence(synthetic_result):
    bg = synthetic_result.tables["background_summary"]
    proxy = bg[bg.background_method == "proxy_stack_3d"].set_index("concentration")
    assert proxy.loc["1:2000", "bg_qc_flags"] != "none"
    assert (proxy["interpretation_note"].str.contains("NOT a measurement of autofluorescence")).all()
    w = " ".join(synthetic_result.warnings)
    assert "RISES as the primary is diluted" in w
    assert "Do not interpret it as autofluorescence" in w


# ------------------------------------------------------------------- z-range modes
def test_manual_z_range_and_bounds():
    s = make_stack(np.ones((10, 8, 8), np.uint16), "1:250", name="a")
    tissue = {"a": np.ones((8, 8), bool)}
    a = AnalysisSettings(depth=DepthSettings(z_range="manual", manual_z_ranges={"a": [2, 7]}))
    ranges, notes = resolve_z_ranges([s], a, tissue)
    assert ranges["a"] == (2, 7) and notes
    bad = AnalysisSettings(depth=DepthSettings(z_range="manual", manual_z_ranges={"a": [2, 12]}))
    with pytest.raises(ValueError, match="outside"):
        resolve_z_ranges([s], bad, tissue)
    with pytest.raises(ValueError, match="no range"):
        resolve_z_ranges([s], AnalysisSettings(depth=DepthSettings(z_range="manual")), tissue)


# --------------------------------------------------------------------- review hash
def _cfg(**analysis):
    return ProjectConfig(series=[SeriesSpec("A", "1:250")], analysis=AnalysisSettings(**analysis))


def test_review_hash_tracks_reviewed_settings_only():
    base = review_settings_sha256(_cfg(), "abc")
    assert base == review_settings_sha256(_cfg(), "abc")
    assert base != review_settings_sha256(_cfg(), "abd")  # different file
    assert base != review_settings_sha256(_cfg(tissue_mask=TissueMaskSettings(method="otsu")), "abc")
    assert base != review_settings_sha256(_cfg(depth=DepthSettings(z_range="auto")), "abc")
    # settings that do not affect masking / z-range do not invalidate a review
    assert base == review_settings_sha256(_cfg(display_percentile=99.5), "abc")


def test_review_status_requires_both_flags_and_matching_hash():
    cfg = _cfg()
    h = review_settings_sha256(cfg, "abc")
    ok, reasons = review_status(cfg, h)
    assert not ok and len(reasons) == 3
    cfg.review.tissue_mask_confirmed = cfg.review.z_range_confirmed = True
    cfg.review.reviewed_settings_sha256 = h
    assert review_status(cfg, h) == (True, [])
    cfg.review.z_range_confirmed = False
    assert not review_status(cfg, h)[0]


# ------------------------------------------------------------- origin separation
def test_data_origin_labels_and_mixing_refused():
    synth = make_concentration_series(shape=(6, 16, 16), include_deconvolved=False)
    assert data_origin(synth) == "synthetic"
    lif = make_stack(np.ones((6, 16, 16), np.uint16), "1:250", name="x")
    object.__setattr__(lif, "source", "curve.lif::x")
    assert data_origin([lif]) == "lif"
    with pytest.raises(ValueError, match="different origin"):
        data_origin(synth + [lif])


def test_every_table_carries_data_origin(synthetic_result):
    for name, df in synthetic_result.tables.items():
        assert (df["data_origin"] == "synthetic").all(), name


def test_demo_output_separated_and_marked(tmp_path):
    with pytest.raises(SystemExit, match="synthetic"):
        main(["demo", "--out", str(tmp_path / "demo_out")])
    out = tmp_path / "synthetic" / "demo"
    assert main(["demo", "--out", str(out)]) == 0
    man = json.loads((out / "run_manifest.json").read_text())
    assert man["synthetic_data"] is True and man["data_origin"] == "synthetic"
    assert "SYNTHETIC DATA" in (out / "summary.md").read_text()


# ------------------------------------------------------------------- plot helpers
def test_element_word():
    assert element_word("stack_all_voxels") == "voxels"
    assert element_word("stack_tissue") == "voxels"
    assert element_word("mip_all_pixels") == "pixels"


@pytest.mark.parametrize("fracs", [[0.50, 0.51, 0.52, 0.53], [0.99, 0.98, 1.0], [0.1, 0.5, 0.9], [0.0, 0.0]])
def test_spread_positions_separates_and_keeps_order(fracs):
    pos = spread_positions(fracs, min_sep=0.07)
    order_in, order_out = np.argsort(fracs, kind="stable"), np.argsort(pos, kind="stable")
    assert list(order_in) == list(order_out)
    sp = np.sort(pos)
    assert np.all(np.diff(sp) >= 0.07 - 1e-9)
    assert sp[-1] <= 1.0 + 1e-9
