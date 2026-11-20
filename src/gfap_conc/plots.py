"""Figures. Every figure states: raw vs deconvolved, display range / transform,
intensity units, depth units and concentration. Synthetic runs are watermarked.

Concentration colours are a single-hue ordinal ramp (darker = more antibody),
validated for lightness order and surface contrast; each concentration also has
its own marker shape and a direct label, so identity never depends on colour alone.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

from .pipeline import AnalysisResult  # noqa: E402
from .projection import shared_display_vmax  # noqa: E402

CONC_COLORS = {"1:250": "#0d366b", "1:500": "#1c5cab", "1:1000": "#3987e5", "1:2000": "#86b6ef"}
CONC_MARKERS = {"1:250": "o", "1:500": "s", "1:1000": "^", "1:2000": "D"}
INK, INK_2, GRID = "#1f1f1e", "#5c5c58", "#e6e6e3"
MASK_BG, MASK_FG = "#2a78d6", "#eb6834"


def _style() -> None:
    plt.rcParams.update({
        "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9, "figure.dpi": 100,
        "savefig.dpi": 200, "savefig.bbox": "tight", "axes.spines.top": False,
        "axes.spines.right": False, "axes.edgecolor": INK_2, "axes.labelcolor": INK,
        "xtick.color": INK_2, "ytick.color": INK_2, "text.color": INK, "axes.grid": True,
        "grid.color": GRID, "grid.linewidth": 0.6, "lines.linewidth": 1.6, "lines.markersize": 5,
        "legend.frameon": False, "figure.facecolor": "white", "axes.facecolor": "white",
    })


def _color(c: str) -> str:
    return CONC_COLORS.get(c, "#7a7a76")


def _marker(c: str) -> str:
    return CONC_MARKERS.get(c, "o")


def _finish(fig, path: Path, synthetic: bool) -> Path:
    if synthetic:
        fig.text(0.5, 0.5, "SYNTHETIC DATA - not experimental results", ha="center", va="center",
                 fontsize=22, color="#d03b3b", alpha=0.22, rotation=20, transform=fig.transFigure, zorder=100)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    plt.close(fig)
    return path


def _scalebar(ax, pixel_size_um: float | None, shape: tuple[int, int]) -> None:
    """Round-number scale bar no longer than ~30% of the image width."""
    if not pixel_size_um:
        return
    ny, nx = shape
    width_um = nx * pixel_size_um
    nice = [v for v in (1, 2, 5, 10, 20, 25, 50, 100, 200, 500) if v <= 0.3 * width_um]
    if not nice:
        return
    length_um = nice[-1]
    n_px = length_um / pixel_size_um
    x0, y0 = nx * 0.95 - n_px, ny * 0.94
    ax.plot([x0, x0 + n_px], [y0, y0], color="white", lw=3, solid_capstyle="butt")
    ax.text(x0 + n_px / 2, y0 - ny * 0.03, f"{length_um:g} µm", color="white", ha="center", va="bottom", fontsize=8)


# ---------------------------------------------------------------- Analysis 1
def mip_grid(result: AnalysisResult, path: Path, scale: str, display_percentile: float, synthetic: bool) -> Path:
    """``scale="shared"``: identical display range for all panels (honest brightness comparison).
    ``scale="per_image"``: each panel 0..its own p99.8 (Rose's display style; hides brightness differences)."""
    _style()
    stacks = result.raw_stacks
    mips = [result.images[s.series_name]["mip"] for s in stacks]
    shared = shared_display_vmax(mips, display_percentile)
    fig, axes = plt.subplots(1, len(stacks), figsize=(3.3 * len(stacks), 3.6), squeeze=False, layout="constrained")
    im = None
    for ax, s, mip in zip(axes[0], stacks, mips, strict=True):
        vmax = shared if scale == "shared" else float(np.percentile(mip, display_percentile))
        im = ax.imshow(mip, cmap="magma", vmin=0, vmax=vmax, interpolation="nearest")
        rng = "" if scale == "shared" else f"\ndisplay 0–{vmax:.0f} (p{display_percentile:g})"
        ax.set_title(f"{s.concentration}  (raw)\n{s.series_name}{rng}", fontsize=8.5)
        ax.set_xticks([]), ax.set_yticks([]), ax.grid(False)
        _scalebar(ax, s.pixel_size_um[1] if s.pixel_size_um else None, mip.shape)
        if scale != "shared":
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    if scale == "shared":
        cb = fig.colorbar(im, ax=axes[0].tolist(), fraction=0.02, pad=0.01, aspect=30)
        cb.set_label(f"{stacks[0].intensity_unit}", fontsize=8)
        sub = f"SHARED display 0–{shared:.0f} for every panel (= largest per-image p{display_percentile:g}); brightness is comparable"
    else:
        sub = (f"PER-IMAGE display 0–p{display_percentile:g} of each MIP ({stacks[0].intensity_unit}); panels look "
               "similarly bright by construction - do not compare brightness here")
    fig.suptitle(f"Raw GFAP max-intensity projections (all z-slices), magma\n{sub}", fontsize=9.5)
    return _finish(fig, path, synthetic)


# ---------------------------------------------------------------- Analysis 2
def depth_profile(result: AnalysisResult, path: Path, metric: str, ylabel: str, synthetic: bool) -> Path | None:
    _style()
    df = result.tables["depth_profile"]
    if metric not in df.columns or df[metric].isna().all():
        return None
    use_um = df["z_um"].notna().all()
    xcol = "z_um" if use_um else "z_index"
    fig, ax = plt.subplots(figsize=(6.2, 3.9))
    for s in result.raw_stacks:
        d = df[df["series_name"] == s.series_name].sort_values("z_index")
        c = s.concentration
        valid = d.dropna(subset=[metric])
        if valid.empty:
            continue
        ax.plot(d[xcol], d[metric], color=_color(c), marker=_marker(c), markevery=max(1, len(d) // 12), label=c)
        last = valid.iloc[-1]
        ax.annotate(c, (last[xcol], last[metric]), xytext=(4, 0), textcoords="offset points",
                    color=INK_2, fontsize=8, va="center")
    ax.set_xlabel("Depth from first acquired slice (µm, nominal)" if use_um else "z-slice index (first acquired = 0)")
    ax.set_ylabel(ylabel)
    region = df["region"].iloc[0]
    ax.set_title(f"Raw GFAP: {ylabel.split(' (')[0].lower()} per z-slice ({region.replace('_', ' ')})", loc="left")
    ax.legend(title="Primary dilution", loc="upper left", bbox_to_anchor=(1.08, 1))
    ax.set_ylim(bottom=0)
    return _finish(fig, path, synthetic)


# ---------------------------------------------------------------- Analysis 3
def background_comparison(result: AnalysisResult, path: Path, synthetic: bool) -> Path | None:
    _style()
    df = result.tables["background_summary"]
    if df.empty:
        return None
    d = df[df["background_method"] == "proxy_stack_3d"]
    if d.empty:
        d = df[df["background_method"] == "manual_roi_stack"]
        if d.empty:
            return None
    method = d["background_method"].iloc[0]
    order = [s.concentration for s in result.raw_stacks if s.concentration in set(d["concentration"])]
    d = d.set_index("concentration").loc[order]
    x = np.arange(len(order))
    panels = [("bg_median", "Background median (raw counts)", "bg_robust_sd"),
              ("signal_to_background", "Foreground median / background median", None)]
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.4))
    for ax, (col, label, err) in zip(axes, panels, strict=True):
        if col not in d or d[col].isna().all():
            ax.set_visible(False)
            continue
        for i, c in enumerate(order):
            ax.errorbar(i, d.loc[c, col], yerr=d.loc[c, err] if err else None, fmt=_marker(c), color=_color(c),
                        markersize=7, capsize=3, lw=1.2)
            ax.annotate(f"{d.loc[c, col]:.3g}", (i, d.loc[c, col]), xytext=(8, 0), textcoords="offset points",
                        fontsize=8, color=INK_2, va="center")
        ax.set_xticks(x, order)
        ax.set_xlim(-0.5, len(order) - 0.5)
        ax.set_xlabel("Primary dilution")
        ax.set_ylabel(label)
        ax.set_ylim(bottom=0)
    axes[0].set_title("Background floor (error bar: robust SD)", loc="left")
    axes[1].set_title("Proxy signal-to-background", loc="left")
    fig.suptitle(f"Raw GFAP background - {method}. BACKGROUND PROXY, not a measured autofluorescence level",
                 fontsize=9.5, y=1.04)
    return _finish(fig, path, synthetic)


def background_overlays(result: AnalysisResult, path: Path, display_percentile: float, synthetic: bool) -> Path | None:
    _style()
    stacks = [s for s in result.raw_stacks if "bg_proxy_mip" in result.images[s.series_name]]
    if not stacks:
        return None
    mips = [result.images[s.series_name]["mip"] for s in stacks]
    vmax = shared_display_vmax(mips, display_percentile)
    fig, axes = plt.subplots(1, len(stacks), figsize=(3.1 * len(stacks), 3.5), squeeze=False)
    for ax, s, mip in zip(axes[0], stacks, mips, strict=True):
        im = result.images[s.series_name]
        g = np.clip(mip / vmax, 0, 1)
        rgb = np.dstack([g, g, g])
        for mask, hexc, alpha in ((im["bg_proxy_mip"], MASK_BG, 0.45), (im["fg_mip"], MASK_FG, 0.55)):
            col = np.array(matplotlib.colors.to_rgb(hexc))
            rgb[mask] = (1 - alpha) * rgb[mask] + alpha * col
        ax.imshow(rgb, interpolation="nearest")
        ax.set_title(f"{s.concentration}  (raw MIP)", fontsize=9)
        ax.set_xticks([]), ax.set_yticks([]), ax.grid(False)
    fig.legend(handles=[Patch(color=MASK_FG, label="GFAP foreground (shared threshold)"),
                        Patch(color=MASK_BG, label="Background proxy region")],
               loc="lower center", ncol=2, bbox_to_anchor=(0.5, -0.06))
    thr = result.params.get("foreground_threshold")
    fig.suptitle(f"Mask QA on raw MIPs (gray: shared display 0–{vmax:.0f} counts; shared threshold {thr:.0f} counts)\n"
                 "2-D MIP masks shown; primary background numbers use slice-wise 3-D masks. Check that faint "
                 "processes are not inside the background region.", fontsize=9.5, y=1.05)
    return _finish(fig, path, synthetic)


# ---------------------------------------------------------------- Analysis 4
def histograms(result: AnalysisResult, path: Path, source: str, zoom_max: float | None, synthetic: bool) -> Path | None:
    _style()
    df = result.tables["intensity_histogram"]
    df = df[df["source"] == source]
    if df.empty:
        return None
    if zoom_max is None:
        mp = result.tables["max_projection_summary"]
        mp = mp[(mp["projection"] == "mip_full_stack") & (mp["region"] == "all_pixels")]
        zoom_max = 1.2 * float(mp["p99_8"].max())
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.7))
    for ax, zoom in zip(axes, (False, True), strict=True):
        for s in result.raw_stacks:
            d = df[df["series_name"] == s.series_name]
            c = s.concentration
            if zoom:
                d = d[d["bin_left"] <= zoom_max]
            ax.step(d["bin_left"], d["density_per_level"], where="post", color=_color(c), label=c, lw=1.4)
        ax.set_yscale("log")
        ax.set_xlabel("Raw intensity (counts)")
        ax.set_ylabel("Fraction of pixels per intensity level (log)")
    w = result.params["histogram_bin_width"]
    axes[0].set_title(f"Full range 0–{result.params['histogram_max_value']} (bin width {w})", loc="left")
    axes[1].set_title(f"Zoom 0–{zoom_max:.0f} (same bins)", loc="left")
    axes[1].legend(title="Primary dilution", loc="upper left", bbox_to_anchor=(1.02, 1))
    label = {"mip_all_pixels": "max-projection pixels", "stack_all_voxels": "all voxels of the z-stack",
             "mip_tissue": "max-projection pixels within tissue mask",
             "stack_tissue": "z-stack voxels within tissue mask"}.get(source, source)
    fig.suptitle(f"Raw GFAP intensity distribution - {label}; identical integer-aligned bins for all concentrations",
                 fontsize=9.5, y=1.03)
    return _finish(fig, path, synthetic)


def survival(result: AnalysisResult, path: Path, source: str, synthetic: bool) -> Path | None:
    _style()
    df = result.tables["intensity_survival"]
    df = df[df["source"] == source]
    if df.empty:
        return None
    fig, ax = plt.subplots(figsize=(6.2, 3.9))
    xmax = 0.0
    for s in result.raw_stacks:
        d = df[df["series_name"] == s.series_name]
        c = s.concentration
        ax.plot(d["threshold"], d["fraction_at_or_above"], color=_color(c), label=c,
                marker=_marker(c), markevery=max(1, len(d) // 20))
        nz = d[d["fraction_at_or_above"] > 1e-4]
        xmax = max(xmax, float(nz["threshold"].max()) if len(nz) else 0)
    ax.set_xlim(0, xmax * 1.05 if xmax else None)
    ax.set_yscale("log")
    ax.set_xlabel("Intensity threshold T (raw counts)")
    ax.set_ylabel("Fraction of pixels ≥ T (log)")
    ax.set_title(f"Raw GFAP survival curves ({source.replace('_', ' ')})", loc="left")
    ax.legend(title="Primary dilution", loc="upper left", bbox_to_anchor=(1.02, 1))
    return _finish(fig, path, synthetic)


def make_all(result: AnalysisResult, fig_dir: Path, display_percentile: float, zoom_max: float | None,
             synthetic: bool = False) -> list[Path]:
    out = [
        mip_grid(result, fig_dir / "raw_mip_grid_shared_scale.png", "shared", display_percentile, synthetic),
        mip_grid(result, fig_dir / "raw_mip_grid_p998_display.png", "per_image", display_percentile, synthetic),
        depth_profile(result, fig_dir / "depth_profile_raw_p998.png", "p99_8",
                      "Per-slice 99.8th percentile (raw counts)", synthetic),
        depth_profile(result, fig_dir / "depth_profile_raw_mean.png", "mean", "Per-slice mean (raw counts)", synthetic),
        depth_profile(result, fig_dir / "depth_profile_background_proxy.png", "bg_proxy_median",
                      "Background-proxy median (raw counts)", synthetic),
        background_comparison(result, fig_dir / "background_comparison.png", synthetic),
        background_overlays(result, fig_dir / "background_ROI_overlays.png", display_percentile, synthetic),
    ]
    for src in ("mip_all_pixels", "stack_all_voxels", "mip_tissue", "stack_tissue"):
        out.append(histograms(result, fig_dir / f"histograms_{src}_shared_bins.png", src, zoom_max, synthetic))
        out.append(survival(result, fig_dir / f"survival_{src}.png", src, synthetic))
    return [p for p in out if p is not None]

