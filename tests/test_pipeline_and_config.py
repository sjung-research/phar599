"""End-to-end pipeline on synthetic data, config validation, and the run CLI on a test LIF."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml
from conftest import make_stack
from lif_writer import LifTestImage, write_lif

from gfap_conc.cli import main
from gfap_conc.config import AnalysisSettings, BackgroundSettings, load_config
from gfap_conc.model import ImageKind
from gfap_conc.pipeline import run_analysis
from gfap_conc.synthetic import make_concentration_series


@pytest.fixture(scope="module")
def result():
    return run_analysis(make_concentration_series(shape=(15, 64, 64)), AnalysisSettings())


def test_deconvolved_excluded_everywhere(result):
    assert len(result.excluded_stacks) == 4
    for name, df in result.tables.items():
        assert "image_kind" in df.columns, name
        assert set(df["image_kind"]) == {"raw"}, name
        assert not df["series_name"].str.contains("_lng").any(), name
    assert any("Excluded" in w for w in result.warnings)


def test_every_row_keeps_identity(result):
    for name, df in result.tables.items():
        for col in ("concentration", "series_name", "channel_index"):
            assert df[col].notna().all(), (name, col)


def test_known_brightness_order_preserved(result):
    mp = result.tables["max_projection_summary"]
    mp = mp[(mp.projection == "mip_full_stack")].set_index("concentration")
    order = ["1:250", "1:500", "1:1000", "1:2000"]
    assert list(mp.loc[order, "mean"]) == sorted(mp["mean"], reverse=True)


def test_mid_depth_dip_detected(result):
    r = result.tables["middle_vs_outer_summary"]
    r = r[r.metric == "p99_8"]
    assert (r["middle_over_outer_mean"] < 0.9).all()


def test_histogram_totals_match_pixel_counts(result):
    h = result.tables["intensity_histogram"]
    for (_, _), g in h.groupby(["series_name", "source"]):
        assert g["count"].sum() == g["n_pixels"].iloc[0]
    assert h.groupby(["series_name", "source"])["bin_left"].apply(tuple).nunique() == 1  # shared bins


def test_background_rows_labelled_as_proxy(result):
    bg = result.tables["background_summary"]
    assert set(bg["background_method"]) == {"proxy_stack_3d", "proxy_mip_2d (inflated by max over z)"}
    s = bg[bg.background_method == "proxy_stack_3d"].set_index("concentration")
    # True synthetic background is 60 counts. The low percentile recovers it everywhere ...
    assert s["bg_p05"].between(40, 70).all()
    assert s.loc["1:250", "bg_median"] == pytest.approx(60, abs=10)
    # ... but with ONE shared threshold, faint processes of the dimmest stain fall below it and
    # are counted as background: a known upward bias of the proxy median for weak stains.
    assert s.loc["1:2000", "bg_median"] > s.loc["1:250", "bg_median"]


def test_missing_background_region_is_warned():
    data = np.full((6, 16, 16), 1000, np.uint16)
    data[:, :, :8] = 10
    a = make_stack(data, "1:250", name="a")
    res = run_analysis([a], AnalysisSettings(background=BackgroundSettings(exclusion_radius_px=20)))
    assert any("EMPTY_REGION" in w and "autofluorescence" in w for w in res.warnings)


def test_duplicate_concentration_not_pooled():
    a = make_stack(np.ones((6, 8, 8), np.uint16), "1:250", name="a")
    b = make_stack(np.ones((6, 8, 8), np.uint16), "1:250", name="b")
    with pytest.raises(ValueError, match="pooled"):
        run_analysis([a, b], AnalysisSettings())


def test_only_deconvolved_refused():
    d = make_stack(np.ones((6, 8, 8), np.uint16), "1:250", kind=ImageKind.DECONVOLVED, bit_depth=16)
    with pytest.raises(ValueError, match="no raw"):
        run_analysis([d], AnalysisSettings())


def test_comparability_warning_for_unequal_z_step():
    a = make_stack(np.ones((6, 8, 8), np.uint16), "1:250", z_step_um=0.5, name="a")
    b = make_stack(np.ones((6, 8, 8), np.uint16), "1:500", z_step_um=1.0, name="b")
    res = run_analysis([a, b], AnalysisSettings())
    assert any("z-step" in w for w in res.warnings)
    proj = set(res.tables["max_projection_summary"]["projection"])
    assert {"mip_matched_depth_first", "mip_matched_depth_center"} <= proj


def test_config_rejects_unknown_keys(tmp_path):
    p = tmp_path / "c.yaml"
    p.write_text(yaml.safe_dump({"analysis": {"background": {"exclusion_radius": 3}}}), encoding="utf-8")
    with pytest.raises(ValueError, match="unknown key"):
        load_config(p)


def test_config_roundtrip(tmp_path):
    p = tmp_path / "c.yaml"
    p.write_text(yaml.safe_dump({
        "input": {"lif_path": "x.lif", "gfap_channel": 0},
        "series": [{"name": "A", "concentration": "1:250"},
                   {"name": "A_lng", "concentration": "1:250", "kind": "deconvolved", "pair_of": "A"}],
        "analysis": {"foreground": {"threshold": 300}},
    }), encoding="utf-8")
    cfg = load_config(p)
    assert cfg.series[1].kind is ImageKind.DECONVOLVED
    assert cfg.analysis.foreground.threshold == 300
    assert cfg.analysis.display_percentile == 99.8


def test_cli_inspect_review_run_on_test_lif(tmp_path, rng):
    images = []
    for i, conc in enumerate(("250", "500", "1000", "2000")):
        base = rng.poisson(50 + 400 / (i + 1), size=(9, 2, 32, 32)).astype(np.uint16)
        images.append(LifTestImage(f"GFAP 1-{conc}", base))
        images.append(LifTestImage(f"GFAP 1-{conc}_lng", base, resolution_bits=16))
    lif = write_lif(tmp_path / "curve.lif", images)
    before = lif.read_bytes()

    assert main(["inspect", str(lif), "--out", str(tmp_path / "inspect")]) == 0
    inv = pd.read_csv(tmp_path / "inspect" / "series_inventory.csv")
    assert len(inv) == 16
    suggestion = yaml.safe_load((tmp_path / "inspect" / "config_suggestion.yaml").read_text(encoding="utf-8"))
    assert suggestion["input"]["gfap_channel"] is None  # never guessed

    cfg = {"input": {"lif_path": str(lif), "lif_sha256": suggestion["input"]["lif_sha256"], "gfap_channel": 0},
           "series": suggestion["series"], "output_dir": str(tmp_path / "results")}
    cfg_path = tmp_path / "cfg.yaml"
    cfg_path.write_text(yaml.safe_dump(cfg), encoding="utf-8")

    # 1. no review -> comparison refused
    with pytest.raises(SystemExit, match="must be reviewed"):
        main(["run", "--config", str(cfg_path), "--out", str(tmp_path / "run0")])

    # 2. explicit preliminary run is allowed but marked everywhere
    assert main(["run", "--config", str(cfg_path), "--out", str(tmp_path / "prelim"), "--allow-unreviewed"]) == 0
    prelim = tmp_path / "prelim"
    assert "PRELIMINARY" in (prelim / "summary.md").read_text(encoding="utf-8")
    man = json.loads((prelim / "run_manifest.json").read_text(encoding="utf-8"))
    assert man["review_status"]["status"].startswith("PRELIMINARY") and man["data_origin"] == "lif"
    assert (prelim / "qc_warnings.txt").read_text(encoding="utf-8").startswith("PRELIMINARY")
    assert "deconvolved series not loaded" in (prelim / "qc_warnings.txt").read_text(encoding="utf-8")

    # 3. review -> record confirmation with hash -> confirmed run
    assert main(["review", "--config", str(cfg_path), "--out", str(tmp_path / "review")]) == 0
    rv = tmp_path / "review"
    for f in ("review_summary.csv", "review_checklist.md", "figures/review_tissue_mask.png",
              "figures/review_z_range.png"):
        assert (rv / f).exists(), f
    assert not (rv / "max_projection_summary.csv").exists()  # review produces no comparison numbers
    h = json.loads((rv / "review_provenance.json").read_text(encoding="utf-8"))["review_settings_sha256"]
    assert h in (rv / "review_checklist.md").read_text(encoding="utf-8")
    cfg["review"] = {"tissue_mask_confirmed": True, "z_range_confirmed": True, "reviewed_settings_sha256": h}
    cfg_path.write_text(yaml.safe_dump(cfg), encoding="utf-8")
    assert main(["run", "--config", str(cfg_path), "--out", str(tmp_path / "run")]) == 0
    out = tmp_path / "run"
    for f in ("max_projection_summary.csv", "depth_profile.csv", "background_summary.csv",
              "threshold_sensitivity.csv", "intensity_histogram.csv", "summary.md", "run_manifest.json",
              "config_used.yaml", "figures/threshold_sensitivity.png"):
        assert (out / f).exists(), f
    assert "PRELIMINARY" not in (out / "summary.md").read_text(encoding="utf-8")
    assert json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))["review_status"]["status"] == "confirmed"
    mp = pd.read_csv(out / "max_projection_summary.csv")
    assert set(mp["image_kind"]) == {"raw"} and len(set(mp["concentration"])) == 4
    assert set(mp["data_origin"]) == {"lif"}

    # 4. changing a reviewed setting invalidates the confirmation
    cfg["analysis"] = {"depth": {"z_range": "manual", "manual_z_ranges": {
        s["name"]: [1, 7] for s in cfg["series"] if s["kind"] == "raw"}}}
    cfg_path.write_text(yaml.safe_dump(cfg), encoding="utf-8")
    with pytest.raises(SystemExit, match="changed since the review"):
        main(["run", "--config", str(cfg_path), "--out", str(tmp_path / "run_changed")])

    # 5. real results may not be written into a 'synthetic' path
    cfg.pop("analysis")
    cfg_path.write_text(yaml.safe_dump(cfg), encoding="utf-8")
    with pytest.raises(SystemExit, match="synthetic"):
        main(["run", "--config", str(cfg_path), "--out", str(tmp_path / "synthetic_mix")])

    assert lif.read_bytes() == before  # source file untouched throughout

    cfg["input"]["lif_sha256"] = "0" * 64
    cfg_path.write_text(yaml.safe_dump(cfg), encoding="utf-8")
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        main(["run", "--config", str(cfg_path), "--out", str(tmp_path / "run2")])


def test_example_config_is_valid():
    cfg = load_config(Path(__file__).resolve().parents[1] / "config" / "analysis.example.yaml")
    assert cfg.review.tissue_mask_confirmed is False
    assert cfg.input.gfap_channel is None


def test_identical_maxima_far_below_ceiling_warned():
    stacks = []
    for i, c in enumerate(("1:250", "1:500", "1:1000")):
        d = np.full((6, 16, 16), 50 + 10 * i, np.uint16)
        d[:, :2, :] = 1420 - i  # common effective ceiling
        stacks.append(make_stack(d, c, bit_depth=16, name=c))
    res = run_analysis(stacks, AnalysisSettings())
    assert any("HARD clipping" in w for w in res.warnings)
    qc = res.tables["qc_stack_summary"]
    assert (qc["n_at_observed_max"] == 6 * 2 * 16).all()


def test_identical_maxima_without_pileup_reports_soft_compression_caveat():
    stacks = []
    for i, c in enumerate(("1:250", "1:500")):
        d = np.full((6, 16, 16), 50, np.uint16)
        d[0, 0, 0] = 1429 - 8 * i  # single brightest voxel, no pile-up
        stacks.append(make_stack(d, c, bit_depth=16, name=c))
    res = run_analysis(stacks, AnalysisSettings())
    w = [x for x in res.warnings if "nearly identical" in x]
    assert w and "no hard clipping" in w[0] and "NOT excluded" in w[0]
