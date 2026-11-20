"""Run Rose's four analyses on a set of GFAP stacks and collect tidy result tables.

The pipeline is data-source agnostic: it receives validated :class:`GFAPStack`
objects (from a .lif via :mod:`gfap_conc.io_lif`, or synthetic stacks in tests) and
returns pandas tables. Deconvolved (``_lng``) stacks may be passed in but are
excluded from every intensity table and listed in the warnings.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from .background import (
    FLAG_MEANING,
    PROXY_NOTE,
    ROI_NOTE,
    background_proxy_mask,
    background_qc_flags,
    background_stats,
    roi_mask,
    signal_metrics,
    threshold_dependence,
    threshold_sensitivity_rows,
)
from .config import AnalysisSettings
from .depth import detect_signal_z_extent, middle_vs_outer, slice_profile
from .histograms import histogram_table, integer_bin_edges, level_counts, quantiles_from_levels, survival_table
from .masks import foreground_mask, pooled_otsu_threshold, tissue_mask_2d
from .model import GFAPStack, ImageKind, sort_by_dilution
from .projection import matched_depth_windows, max_projection, projection_rows, window_projection
from .stats import summarize


@dataclass
class AnalysisResult:
    tables: dict[str, pd.DataFrame]
    warnings: list[str]
    raw_stacks: list[GFAPStack]
    excluded_stacks: list[GFAPStack]
    images: dict[str, dict[str, np.ndarray]] = field(default_factory=dict)
    params: dict = field(default_factory=dict)


def comparability_warnings(raw: list[GFAPStack]) -> list[str]:
    """Differences between raw stacks that limit absolute-intensity comparison."""
    out = []

    def check(label: str, values: dict[str, object]) -> None:
        if len({repr(v) for v in values.values()}) > 1:
            detail = ", ".join(f"{k}={v}" for k, v in values.items())
            out.append(f"Raw stacks differ in {label}: {detail}")

    check("bit depth", {s.concentration: s.bit_depth for s in raw})
    check("dtype", {s.concentration: str(s.data.dtype) for s in raw})
    check("XY size", {s.concentration: s.data.shape[1:] for s in raw})
    check("number of z-slices", {s.concentration: s.n_z for s in raw})
    check("z-step (um)", {s.concentration: s.z_step_um for s in raw})
    check("XY pixel size (um)", {s.concentration: s.pixel_size_um for s in raw})
    check("GFAP channel index", {s.concentration: s.channel_index for s in raw})
    check("intensity unit", {s.concentration: s.intensity_unit for s in raw})
    if any(s.z_step_um is None for s in raw):
        out.append("At least one raw stack has no z-step; physical depth is unavailable.")
    return out


def _effective_mask(mask: np.ndarray) -> np.ndarray | None:
    return None if mask.all() else mask


def data_origin(stacks: list[GFAPStack]) -> str:
    """``synthetic`` | ``lif`` | ``unspecified``; mixing origins in one analysis is refused."""
    def one(s: GFAPStack) -> str:
        if s.source.startswith("synthetic"):
            return "synthetic"
        return "lif" if "::" in s.source else "unspecified"

    origins = {one(s) for s in stacks}
    if len(origins) > 1:
        raise ValueError(f"stacks of different origin {sorted(origins)} cannot be analysed together")
    return origins.pop()


def compute_tissue_masks(
    raw: list[GFAPStack],
    settings: AnalysisSettings,
    reference_images: dict[str, np.ndarray] | None = None,
) -> dict[str, np.ndarray]:
    """2-D tissue mask per raw stack, exactly as used by the analysis (and by `review`)."""
    t = settings.tissue_mask
    out = {}
    for s in raw:
        if t.reference == "gfap_mip":
            ref = max_projection(s.data)
        else:
            if not reference_images or s.series_name not in reference_images:
                raise ValueError(f"tissue mask needs a reference-channel image for {s.series_name!r}")
            ref = reference_images[s.series_name]
        tm = tissue_mask_2d(ref, t.method, t.threshold, t.sigma_px, t.min_object_px)
        if not tm.any():
            raise ValueError(f"tissue mask for {s.series_name!r} is empty; check tissue_mask settings")
        out[s.series_name] = tm
    return out


def auto_z_range(stack: GFAPStack, tissue_mask: np.ndarray, settings: AnalysisSettings) -> tuple[int, int]:
    prof = slice_profile(stack, _effective_mask(tissue_mask))
    d = settings.depth
    return detect_signal_z_extent(prof[d.auto_metric].to_numpy(), d.auto_rel_level)


def resolve_z_ranges(
    raw: list[GFAPStack], settings: AnalysisSettings, tissue: dict[str, np.ndarray]
) -> tuple[dict[str, tuple[int, int]], list[str]]:
    """Inclusive analysed slice range per stack for z_range = full | auto | manual."""
    d = settings.depth
    ranges, notes = {}, []
    for s in raw:
        if d.z_range == "auto":
            zr = auto_z_range(s, tissue[s.series_name], settings)
        elif d.z_range == "manual":
            if s.series_name not in d.manual_z_ranges:
                raise ValueError(f"z_range is manual but no range is configured for {s.series_name!r}")
            z0, z1 = (int(v) for v in d.manual_z_ranges[s.series_name])
            if not 0 <= z0 <= z1 <= s.n_z - 1:
                raise ValueError(f"{s.series_name}: manual z-range [{z0}, {z1}] outside 0-{s.n_z - 1}")
            zr = (z0, z1)
        else:
            zr = (0, s.n_z - 1)
        if zr != (0, s.n_z - 1):
            notes.append(f"{s.series_name}: {d.z_range} z-range keeps slices {zr[0]}-{zr[1]} of 0-{s.n_z - 1}")
        ranges[s.series_name] = zr
    return ranges, notes


def run_analysis(
    stacks: list[GFAPStack],
    settings: AnalysisSettings,
    reference_images: dict[str, np.ndarray] | None = None,
) -> AnalysisResult:
    raw = sort_by_dilution([s for s in stacks if s.kind is ImageKind.RAW])
    excluded = sort_by_dilution([s for s in stacks if s.kind is not ImageKind.RAW])
    if not raw:
        raise ValueError("no raw stacks supplied; deconvolved (_lng) images cannot answer the concentration question")
    dup = [c for c, n in Counter(s.concentration for s in raw).items() if n > 1]
    if dup:
        raise ValueError(
            f"more than one raw stack for concentration(s) {dup}. Replicates/fields need an explicit "
            "aggregation design (not implemented yet) - they must not be pooled silently."
        )

    warnings = comparability_warnings(raw)
    if excluded:
        warnings.append(
            "Excluded from all intensity analyses (deconvolved, full-range normalised): "
            + ", ".join(s.series_name for s in excluded)
        )
    expected = {"1:250", "1:500", "1:1000", "1:2000"}
    missing = expected - {s.concentration for s in raw}
    if missing:
        warnings.append(f"No raw stack for expected concentration(s): {sorted(missing, key=lambda c: int(c[2:]))}")

    a = settings
    origin = data_origin(raw + excluded)

    # ---------------------------------------------------------------- masks
    tissue = compute_tissue_masks(raw, a, reference_images)
    images: dict[str, dict[str, np.ndarray]] = {
        s.series_name: {"mip": max_projection(s.data), "tissue_mask": tissue[s.series_name]} for s in raw
    }

    thr_spec = a.foreground.threshold
    if thr_spec is None:
        fg_threshold, fg_source = None, "disabled"
        warnings.append("Foreground threshold disabled: foreground fraction, background proxy and SBR not computed.")
    elif isinstance(thr_spec, str):
        arrays = [s.data[:, tissue[s.series_name]] for s in raw]
        fg_threshold = pooled_otsu_threshold(arrays, a.foreground.sample_per_stack, a.foreground.seed)
        fg_source = "pooled_otsu (equal pixel sample per raw stack, shared by all concentrations)"
    else:
        fg_threshold, fg_source = float(thr_spec), "fixed in config (raw counts)"

    # ------------------------------------------------------------ z ranges
    z_ranges, z_notes = resolve_z_ranges(raw, a, tissue)
    warnings += z_notes

    # ------------------------------------------- Analysis 1: projections
    proj_rows: list[dict] = []
    windows = matched_depth_windows(raw, z_ranges)
    # matched-depth rows are only informative if some stack is actually cropped
    emit_matched = any(
        w != (0, s.n_z) for s in raw for w in windows[s.series_name].values()
    )
    for s in raw:
        tm = tissue[s.series_name]
        mip = images[s.series_name]["mip"]
        proj_rows += projection_rows(s, mip, "mip_full_stack", (0, s.n_z), tm, a.display_percentile)
        z0, z1 = z_ranges[s.series_name]
        if (z0, z1) != (0, s.n_z - 1):
            proj_rows += projection_rows(
                s, window_projection(s.data, z0, z1 + 1), "mip_analysis_range", (z0, z1 + 1), tm, a.display_percentile
            )
        for place, (w0, w1) in windows[s.series_name].items() if emit_matched else ():
            proj_rows += projection_rows(
                s, window_projection(s.data, w0, w1), f"mip_matched_depth_{place}", (w0, w1), tm, a.display_percentile
            )

    # ------------------------------------------- Analysis 2: depth profile
    profiles, retention_rows = [], []
    for s in raw:
        tm_eff = _effective_mask(tissue[s.series_name])
        prof = slice_profile(s, tm_eff, fg_threshold, z_ranges[s.series_name])
        profiles.append(prof)
        for metric in a.depth.retention_metrics:
            try:
                r = middle_vs_outer(prof, metric, z_ranges[s.series_name], a.depth.window_fraction)
            except ValueError as exc:
                warnings.append(f"{s.series_name}: middle-vs-outer not computed for {metric}: {exc}")
                continue
            retention_rows.append({**s.identity(), "region": prof["region"].iloc[0], **r})

    # ------------------------------------------- Analysis 3: background
    bg_rows, sens_rows = [], []
    for s, prof in zip(raw, profiles, strict=True):
        tm = tissue[s.series_name]
        tm_eff = _effective_mask(tm)
        z0, z1 = z_ranges[s.series_name]
        data = s.data[z0 : z1 + 1]
        tissue3d = np.broadcast_to(tm, data.shape)
        n_tissue_vox = int(tm.sum()) * data.shape[0]
        common = {**s.identity(), "z_start": z0, "z_stop_exclusive": z1 + 1}
        if fg_threshold is not None:
            r = a.background.exclusion_radius_px
            fg3d = foreground_mask(data, fg_threshold) & tissue3d
            bg3d = background_proxy_mask(fg3d, tm_eff, r)
            bg = background_stats(data, bg3d, n_tissue_vox)

            # threshold-sensitivity diagnostic (shared multiples of the one shared threshold)
            sens = threshold_sensitivity_rows(
                data, tm_eff, [(f, f * fg_threshold) for f in sorted(a.foreground.sensitivity_factors)], r
            )
            sens_rows += [{**common, "exclusion_radius_px": r, **row} for row in sens]
            dep = threshold_dependence(sens)
            flags = background_qc_flags(bg, dep, a.background.qc_max_fraction_of_tissue,
                                        a.background.qc_max_tail_asymmetry, a.background.qc_max_threshold_dependence)
            for f in flags:
                warnings.append(f"{s.series_name} [{s.concentration}] background proxy {f}: {FLAG_MEANING[f]}. "
                                "Do not read this background value as autofluorescence.")
            bg_rows.append(
                {**common, "background_method": "proxy_stack_3d", "foreground_threshold": fg_threshold,
                 "exclusion_radius_px": r, **bg, "bg_threshold_dependence": dep,
                 "bg_qc_flags": ";".join(flags) or "none", "interpretation_note": PROXY_NOTE,
                 **signal_metrics(data, fg3d, bg)}
            )
            mip = images[s.series_name]["mip"]
            fg2d = foreground_mask(mip, fg_threshold) & tm
            bg2d = background_proxy_mask(fg2d, tm_eff, r)
            bgm = background_stats(mip, bg2d, int(tm.sum()))
            bg_rows.append(
                {**common, "background_method": "proxy_mip_2d (inflated by max over z)",
                 "foreground_threshold": fg_threshold, "exclusion_radius_px": r, **bgm,
                 "bg_qc_flags": "visual QA only", "interpretation_note": PROXY_NOTE,
                 **signal_metrics(mip, fg2d, bgm)}
            )
            images[s.series_name]["bg_proxy_mip"] = bg2d
            images[s.series_name]["fg_mip"] = fg2d
            # per-slice background floor, added to the depth profile
            med = np.full(s.n_z, np.nan)
            cnt = np.zeros(s.n_z, dtype=np.int64)
            for i in range(data.shape[0]):
                v = data[i][bg3d[i]]
                cnt[z0 + i] = v.size
                if v.size:
                    med[z0 + i] = float(np.median(v))
            prof["bg_proxy_median"] = med
            prof["bg_proxy_n_pixels"] = cnt
        rois = a.background.manual_rois.get(s.series_name)
        if rois:
            m2 = roi_mask(tm.shape, rois)
            m3 = np.broadcast_to(m2, data.shape)
            bg = background_stats(data, m3, n_tissue_vox)
            row = {**common, "background_method": "manual_roi_stack", "manual_rois": str(rois), **bg,
                   "interpretation_note": ROI_NOTE}
            if fg_threshold is not None:
                row.update(signal_metrics(data, foreground_mask(data, fg_threshold) & tissue3d, bg))
            bg_rows.append(row)

    proxy = [r for r in bg_rows if r["background_method"] == "proxy_stack_3d" and np.isfinite(r["bg_median"])]
    if len(proxy) >= 3:
        meds = [r["bg_median"] for r in sorted(proxy, key=lambda r: r["dilution_factor"])]
        # never decreasing with dilution and >5% higher at the weakest stain than at the strongest
        if all(b >= a_ for a_, b in zip(meds, meds[1:], strict=False)) and meds[-1] > 1.05 * meds[0]:
            warnings.append(
                "Background proxy RISES as the primary is diluted (" + ", ".join(
                    f"{r['concentration']}={r['bg_median']:.4g}" for r in sorted(proxy, key=lambda r: r["dilution_factor"]))
                + "). Less antibody should not create more background; this pattern is expected when faint GFAP "
                "of weaker stains falls below the shared threshold. Do not interpret it as autofluorescence."
            )

    # ------------------------------------------- Analysis 4: histograms
    max_value = max(s.saturation_value for s in raw)
    edges = integer_bin_edges(max_value, a.histogram.n_bins)
    hist_parts, surv_parts, quant_rows = [], [], []
    for s in raw:
        tm = tissue[s.series_name]
        z0, z1 = z_ranges[s.series_name]
        mip = images[s.series_name]["mip"]
        data = s.data[z0 : z1 + 1]
        sources: list[tuple[str, np.ndarray, np.ndarray | None]] = [("mip_all_pixels", mip, None), ("stack_all_voxels", data, None)]
        if not tm.all():
            sources += [("mip_tissue", mip, tm), ("stack_tissue", data, np.broadcast_to(tm, data.shape))]
        for source, arr, mask in sources:
            levels = level_counts(arr, max_value, mask)
            expected_n = arr.size if mask is None else int(np.count_nonzero(mask))
            if int(levels.sum()) != expected_n:
                raise AssertionError(f"histogram total {levels.sum()} != pixel count {expected_n}")
            ident = {**s.identity(), "source": source, "n_pixels": expected_n}
            hist_parts.append(histogram_table(levels, edges, max_value).assign(**ident))
            surv_parts.append(survival_table(levels, edges[:-1]).assign(**ident))
            quant_rows.append({**ident, **quantiles_from_levels(levels)})

    # ------------------------------------------- QC
    qc_rows = []
    for s in raw:
        st = summarize(s.data, s.saturation_value)
        per_slice_sat = [np.count_nonzero(p >= s.saturation_value) / p.size for p in s.data]
        z0, z1 = z_ranges[s.series_name]
        qc_rows.append(
            {
                **s.identity(),
                "source": s.source,
                "shape_zyx": "x".join(map(str, s.data.shape)),
                "dtype": str(s.data.dtype),
                "z_step_um": s.z_step_um,
                "pixel_size_y_um": s.pixel_size_um[0] if s.pixel_size_um else np.nan,
                "pixel_size_x_um": s.pixel_size_um[1] if s.pixel_size_um else np.nan,
                "saturation_value": s.saturation_value,
                "stack_min": st["min"],
                "stack_max": st["max"],
                "stack_mean": st["mean"],
                "stack_frac_zero": st["frac_zero"],
                "stack_frac_saturated": st["frac_saturated"],
                "max_slice_frac_saturated": float(max(per_slice_sat)),
                "z_analysis_first": z0,
                "z_analysis_last": z1,
                "tissue_mask_method": a.tissue_mask.method,
                "tissue_fraction_of_field": float(tissue[s.series_name].mean()),
                "foreground_threshold": fg_threshold,
            }
        )
        if st["frac_saturated"] > 0.001:
            warnings.append(f"{s.series_name}: {st['frac_saturated']:.3%} of voxels are saturated")
        if st["frac_zero"] > 0.01:
            warnings.append(f"{s.series_name}: {st['frac_zero']:.2%} of voxels are exactly 0 (offset clipping?)")

    tables = {
        "qc_stack_summary": pd.DataFrame(qc_rows),
        "max_projection_summary": pd.DataFrame(proj_rows),
        "depth_profile": pd.concat(profiles, ignore_index=True),
        "middle_vs_outer_summary": pd.DataFrame(retention_rows),
        "background_summary": pd.DataFrame(bg_rows),
        "intensity_histogram": pd.concat(hist_parts, ignore_index=True),
        "intensity_survival": pd.concat(surv_parts, ignore_index=True),
        "histogram_quantiles": pd.DataFrame(quant_rows),
        "threshold_sensitivity": pd.DataFrame(sens_rows),
    }
    for df in tables.values():
        df.insert(0, "data_origin", origin)
    params = {
        "data_origin": origin,
        "foreground_threshold": fg_threshold,
        "foreground_threshold_source": fg_source,
        "histogram_max_value": int(max_value),
        "histogram_bin_width": int(edges[1] - edges[0]),
        "z_ranges": {k: list(v) for k, v in z_ranges.items()},
        "matched_depth_windows": {k: {p: list(w) for p, w in v.items()} for k, v in windows.items()},
    }
    return AnalysisResult(tables, warnings, raw, excluded, images, params)


def write_tables(result: AnalysisResult, out_dir: Path) -> list[Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for name, df in result.tables.items():
        p = out_dir / f"{name}.csv"
        df.to_csv(p, index=False)
        paths.append(p)
    p = out_dir / "qc_warnings.txt"
    p.write_text("\n".join(result.warnings) + ("\n" if result.warnings else ""), encoding="utf-8")
    paths.append(p)
    return paths
