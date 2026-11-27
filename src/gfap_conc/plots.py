"""Figures. Every figure states: raw vs deconvolved, display range / transform,
intensity units, depth units and concentration.

Watermarks: synthetic runs carry "SYNTHETIC DATA"; real runs whose tissue mask /
z-range have not been confirmed carry "PRELIMINARY". Reviewed real runs carry none.

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

from .model import GFAPStack  # noqa: E402
from .pipeline import AnalysisResult  # noqa: E402
from .projection import max_projection, shared_display_vmax  # noqa: E402

CONC_COLORS = {"1:250": "#0d366b", "1:500": "#1c5cab", "1:1000": "#3987e5", "1:2000": "#86b6ef"}
CONC_MARKERS = {"1:250": "o", "1:500": "s", "1:1000": "^", "1:2000": "D"}
INK, INK_2, GRID = "#1f1f1e", "#5c5c58", "#e6e6e3"
MASK_BG, MASK_FG = "#2a78d6", "#eb6834"
FLAG_INK = "#d03b3b"

SYNTHETIC_WATERMARK = "SYNTHETIC DATA - not experimental results"
PRELIMINARY_WATERMARK = "PRELIMINARY - tissue mask / z-range not reviewed"


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


def element_word(source: str) -> str:
    """'voxels' for whole-stack data, 'pixels' for projections."""
    return "voxels" if source.startswith("stack") else "pixels"


def _finish(fig, path: Path, watermark: str | None) -> Path:
    if watermark:
        fig.text(0.5, 0.5, watermark, ha="center", va="center", fontsize=20, color=FLAG_INK, alpha=0.22,
                 rotation=20, transform=fig.transFigure, zorder=100)
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


def spread_positions(fracs: list[float], min_sep: float = 0.07, lo: float = 0.0, hi: float = 1.0) -> list[float]:
    """Nudge label positions (axes fractions) apart by at least ``min_sep``, keeping order
    and staying within [lo, hi] where possible. Returns positions in the input order."""
    if not fracs:
        return []
    order = np.argsort(fracs)
    pos = np.asarray(fracs, float)[order].copy()
    for i in range(1, len(pos)):
        pos[i] = max(pos[i], pos[i - 1] + min_sep)
    if pos[-1] > hi:
        pos[-1] = hi
        for i in range(len(pos) - 2, -1, -1):
            pos[i] = min(pos[i], pos[i + 1] - min_sep)
    if pos[0] < lo:  # too many labels for the space: spread evenly
        pos = np.linspace(lo, hi, len(pos))
    out = np.empty_like(pos)
    out[order] = pos
    return out.tolist()


def _end_labels(ax, items: list[tuple[float, float, str]]) -> None:
    """Direct labels just right of the plot area, de-overlapped, with leader lines."""
    if not items:
        return
    y0, y1 = ax.get_ylim()
    span = (y1 - y0) or 1.0
    fr = [(y - y0) / span for _, y, _ in items]
    pos = spread_positions(fr)
    for (x, y, text), f in zip(items, pos, strict=True):
        ax.annotate(text, xy=(x, y), xycoords="data", xytext=(1.02, f), textcoords="axes fraction",
                    va="center", fontsize=8, color=INK_2,
                    arrowprops={"arrowstyle": "-", "color": "#b4b4b0", "lw": 0.7, "shrinkA": 1, "shrinkB": 2})


# ---------------------------------------------------------------- Analysis 1
def mip_grid(result: AnalysisResult, path: Path, scale: str, display_percentile: float,
             watermark: str | None) -> Path:
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
    return _finish(fig, path, watermark)


# ---------------------------------------------------------------- Analysis 2
def depth_profile(result: AnalysisResult, path: Path, metric: str, ylabel: str, watermark: str | None,
                  note: str | None = None) -> Path | None:
    _style()
    df = result.tables["depth_profile"]
    if metric not in df.columns or df[metric].isna().all():
        return None
    use_um = df["z_um"].notna().all()
    xcol = "z_um" if use_um else "z_index"
    fig, ax = plt.subplots(figsize=(6.4, 3.9))
    ends = []
    for s in result.raw_stacks:
        d = df[df["series_name"] == s.series_name].sort_values("z_index")
        c = s.concentration
        valid = d.dropna(subset=[metric])
        if valid.empty:
            continue
        ax.plot(d[xcol], d[metric], color=_color(c), marker=_marker(c), markevery=max(1, len(d) // 12), label=c)
        last = valid.iloc[-1]
        ends.append((float(last[xcol]), float(last[metric]), c))
    ax.set_xlabel("Depth from first acquired slice (µm, nominal)" if use_um else "z-slice index (first acquired = 0)")
    ax.set_ylabel(ylabel)
    ax.set_ylim(bottom=0)
    _end_labels(ax, ends)
    region = df["region"].iloc[0]
    ax.set_title(f"Raw GFAP: {ylabel.split(' (')[0].lower()} per z-slice ({region.replace('_', ' ')})", loc="left")
    ax.legend(title="Primary dilution", loc="upper left", bbox_to_anchor=(1.16, 1))
    if note:
        fig.text(0.01, -0.04, note, fontsize=7.5, color=INK_2, ha="left", va="top", wrap=True)
    return _finish(fig, path, watermark)


# ---------------------------------------------------------------- Analysis 3
ERRORBAR_NOTE = ("Error bars: robust SD (1.4826 × MAD) of background-proxy voxels WITHIN one image - spatial and "
                 "noise spread inside a single field. Not a confidence interval and not biological-replicate "
                 "uncertainty (n = 1 field per concentration).")


def background_comparison(result: AnalysisResult, path: Path, watermark: str | None) -> Path | None:
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
    panels = [("bg_median", "Background-proxy median (raw counts)", "bg_robust_sd"),
              ("signal_to_background", "Foreground median / background median", None)]
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.5))
    any_flag = False
    for ax, (col, label, err) in zip(axes, panels, strict=True):
        if col not in d or d[col].isna().all():
            ax.set_visible(False)
            continue
        for i, c in enumerate(order):
            ax.errorbar(i, d.loc[c, col], yerr=d.loc[c, err] if err else None, fmt=_marker(c), color=_color(c),
                        markersize=7, capsize=3, lw=1.2)
            ax.annotate(f"{d.loc[c, col]:.3g}", (i, d.loc[c, col]), xytext=(8, 0), textcoords="offset points",
                        fontsize=8, color=INK_2, va="center")
            flags = str(d.loc[c].get("bg_qc_flags", "none"))
            if flags not in ("none", "nan"):
                any_flag = True
                ax.annotate("QC", (i, 0), xytext=(0, 4), textcoords="offset points", ha="center", va="bottom",
                            fontsize=7.5, color=FLAG_INK, fontweight="bold")
        ax.set_xticks(x, order)
        ax.set_xlim(-0.5, len(order) - 0.5)
        ax.set_xlabel("Primary dilution")
        ax.set_ylabel(label)
        ax.set_ylim(bottom=0)
    axes[0].set_title("Background floor: median ± robust within-image SD", loc="left")
    axes[1].set_title("Proxy signal-to-background", loc="left")
    fig.suptitle(f"Raw GFAP background - {method}. BACKGROUND PROXY, not a measured autofluorescence level",
                 fontsize=9.5, y=1.04)
    note = ERRORBAR_NOTE
    if any_flag:
        note += ("\n'QC' marks concentrations whose background proxy is flagged (bg_qc_flags in "
                 "background_summary.csv; see qc_warnings.txt). Do not interpret flagged values as autofluorescence.")
    fig.text(0.01, -0.03, note, fontsize=7.5, color=INK_2, ha="left", va="top", wrap=True)
    return _finish(fig, path, watermark)


def background_overlays(result: AnalysisResult, path: Path, display_percentile: float,
                        watermark: str | None) -> Path | None:
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
    return _finish(fig, path, watermark)


def threshold_sensitivity(result: AnalysisResult, path: Path, watermark: str | None) -> Path | None:
    _style()
    df = result.tables.get("threshold_sensitivity")
    if df is None or df.empty:
        return None
    panels = [("bg_median", "Background-proxy median (raw counts)"),
              ("bg_fraction_of_reference", "Proxy fraction of tissue voxels"),
              ("signal_to_background", "Proxy signal-to-background")]
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 3.6), layout="constrained")
    thr = result.params.get("foreground_threshold")
    for ax, (col, label) in zip(axes, panels, strict=True):
        for s in result.raw_stacks:
            d = df[df["series_name"] == s.series_name].sort_values("foreground_threshold")
            c = s.concentration
            ax.plot(d["foreground_threshold"], d[col], color=_color(c), marker=_marker(c), label=c)
        if thr:
            ax.axvline(thr, color=INK_2, lw=0.9, ls="--")
        ax.set_xlabel("Shared foreground threshold T (raw counts)")
        ax.set_ylabel(label)
        ax.set_ylim(bottom=0)
    axes[0].set_title("Stable = proxy is noise-like; rising = contains faint GFAP", loc="left", fontsize=9)
    axes[1].set_title("Share of tissue classed as background", loc="left", fontsize=9)
    axes[2].set_title("How much SBR depends on T", loc="left", fontsize=9)
    axes[2].legend(title="Primary dilution", loc="upper left", bbox_to_anchor=(1.02, 1))
    factors = ", ".join(f"{f:g}" for f in sorted(df["threshold_factor"].unique()))
    fig.suptitle(f"Threshold-sensitivity DIAGNOSTIC (raw stacks, proxy background). Dashed line: threshold used "
                 f"({thr:.0f} counts); others = {factors} × that value, applied to all concentrations alike",
                 fontsize=9.5)
    return _finish(fig, path, watermark)


# ---------------------------------------------------------------- Analysis 4
def histograms(result: AnalysisResult, path: Path, source: str, zoom_max: float | None,
               watermark: str | None) -> Path | None:
    _style()
    df = result.tables["intensity_histogram"]
    df = df[df["source"] == source]
    if df.empty:
        return None
    if zoom_max is None:
        mp = result.tables["max_projection_summary"]
        mp = mp[(mp["projection"] == "mip_full_stack") & (mp["region"] == "all_pixels")]
        zoom_max = 1.2 * float(mp["p99_8"].max())
    word = element_word(source)
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
        ax.set_ylabel(f"Fraction of {word} per intensity level (log)")
    w = result.params["histogram_bin_width"]
    axes[0].set_title(f"0–{result.params['histogram_max_value']} ({result.params.get('histogram_range_basis', '')}; "
                      f"bin width {w})", loc="left", fontsize=8.5)
    axes[1].set_title(f"Zoom 0–{zoom_max:.0f} (same bins)", loc="left")
    axes[1].legend(title="Primary dilution", loc="upper left", bbox_to_anchor=(1.02, 1))
    label = {"mip_all_pixels": "max-projection pixels", "stack_all_voxels": "all voxels of the z-stack",
             "mip_tissue": "max-projection pixels within tissue mask",
             "stack_tissue": "z-stack voxels within tissue mask"}.get(source, source)
    fig.suptitle(f"Raw GFAP intensity distribution - {label}; identical integer-aligned bins for all concentrations",
                 fontsize=9.5, y=1.03)
    return _finish(fig, path, watermark)


def survival(result: AnalysisResult, path: Path, source: str, watermark: str | None) -> Path | None:
    _style()
    df = result.tables["intensity_survival"]
    df = df[df["source"] == source]
    if df.empty:
        return None
    word = element_word(source)
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
    ax.set_ylabel(f"Fraction of {word} ≥ T (log)")
    ax.set_title(f"Raw GFAP survival curves ({source.replace('_', ' ')})", loc="left")
    ax.legend(title="Primary dilution", loc="upper left", bbox_to_anchor=(1.02, 1))
    return _finish(fig, path, watermark)


# ---------------------------------------------------------------- Rose-style figures
def _short_unit(unit: str) -> str:
    return "Photon counts" if "photon" in unit.lower() else "Raw counts"


def mip_with_histograms(result: AnalysisResult, path: Path, display_percentile: float,
                        watermark: str | None) -> Path | None:
    """Rose's first-staining layout: each MIP (per-image display 0..p99.8) above the histogram
    of that MIP's pixels (raw pixel counts, log y) with the p99.8 marked in red.

    Unlike Rose's figure, all histograms share bins, x-range and y-range, so columns are
    directly comparable. Bin counts come from intensity_histogram.csv (source mip_all_pixels).
    """
    _style()
    hist = result.tables["intensity_histogram"]
    hist = hist[hist["source"] == "mip_all_pixels"]
    mp = result.tables["max_projection_summary"]
    mp = mp[(mp["projection"] == "mip_full_stack") & (mp["region"] == "all_pixels")].set_index("series_name")
    stacks = [s for s in result.raw_stacks if s.series_name in set(hist["series_name"])]
    if not stacks:
        return None
    n = len(stacks)
    fig = plt.figure(figsize=(3.7 * n, 5.9))
    gs = fig.add_gridspec(2, n, height_ratios=[1.55, 1], hspace=0.42, wspace=0.28, top=0.86, bottom=0.08)
    hax0 = None
    for j, s in enumerate(stacks):
        mip = result.images[s.series_name]["mip"]
        vmax = float(mp.loc[s.series_name, "display_vmax_within_image"])
        ax = fig.add_subplot(gs[0, j])
        im = ax.imshow(mip, cmap="magma", vmin=0, vmax=vmax, interpolation="nearest")
        ax.set_xticks([]), ax.set_yticks([]), ax.grid(False)
        ax.set_title(f"{s.concentration}  ({s.series_name})\ndisplay 0–{vmax:.0f} counts, "
                     f"{display_percentile:g}th pct", fontsize=8.5)
        _scalebar(ax, s.pixel_size_um[1] if s.pixel_size_um else None, mip.shape)
        cb = fig.colorbar(im, ax=ax, orientation="horizontal", fraction=0.05, pad=0.03, aspect=30)
        cb.set_label(_short_unit(s.intensity_unit), fontsize=7.5)
        cb.ax.tick_params(labelsize=7)

        hax = fig.add_subplot(gs[1, j], sharex=hax0, sharey=hax0)
        hax0 = hax0 or hax
        d = hist[hist["series_name"] == s.series_name]
        d = d[d["count"] > 0]
        hax.bar(d["bin_left"], d["count"], width=d["bin_right_exclusive"] - d["bin_left"], align="edge",
                color="#555553", linewidth=0)
        hax.set_yscale("log")
        hax.axvline(vmax, color=FLAG_INK, lw=1.2)
        x_hi = result.params["histogram_max_value"] + result.params["histogram_bin_width"]
        right = vmax < 0.4 * x_hi  # near the left edge: put the label right of the line
        hax.annotate(f"{display_percentile:g}th = {vmax:.0f}", xy=(vmax, 0.97), xycoords=hax.get_xaxis_transform(),
                     xytext=(4 if right else -4, 0), textcoords="offset points", ha="left" if right else "right",
                     va="top", color=FLAG_INK, fontsize=8)
        hax.set_xlabel(f"Max-projection pixel intensity ({_short_unit(s.intensity_unit).lower()})", fontsize=8)
        if j == 0:
            hax.set_ylabel("Pixels (log)")
        hax.grid(False)
    hax0.set_xlim(0, result.params["histogram_max_value"] + result.params["histogram_bin_width"])
    w = result.params["histogram_bin_width"]
    fig.suptitle("Raw GFAP max projections with their pixel-intensity histograms (layout of Rose's first-staining "
                 f"figure)\nTop: PER-IMAGE display 0–p{display_percentile:g} (brightness not comparable). Bottom: "
                 f"MIP pixel counts, shared {w}-count bins and axes; red line = that image's p{display_percentile:g}",
                 fontsize=9.5)
    return _finish(fig, path, watermark)


def depth_slice_gallery(result: AnalysisResult, path: Path, scale: str, display_percentile: float,
                        watermark: str | None) -> Path | None:
    """Rose's depth layout: one row per stack, XY slices at 0/25/50/75/90% depth and the
    per-slice p99.8 profile with the shown slices marked.

    ``scale="per_image"`` (Rose): all slices of a row share that stack's own display range,
    0..p99.8 of all its voxels. ``scale="shared"``: one range for every row (the largest of
    those), so rows are comparable.
    """
    from .depth import GALLERY_FRACTIONS, gallery_slice_indices

    _style()
    stacks = result.raw_stacks
    prof = result.tables["depth_profile"]
    if not stacks or "p99_8" not in prof:
        return None
    vmaxes = {s.series_name: float(np.percentile(s.data, display_percentile)) for s in stacks}
    shared = max(vmaxes.values())
    n, k = len(stacks), len(GALLERY_FRACTIONS)
    s0 = stacks[0]
    aspect = s0.data.shape[1] / s0.data.shape[2]
    fig = plt.figure(figsize=(2.3 * k + 4.6, 0.25 + n * (2.3 * aspect + 0.75)))
    gs = fig.add_gridspec(n, k + 3, width_ratios=[1] * k + [0.06, 0.42, 1.7], wspace=0.08, hspace=0.45,
                          top=1 - 0.9 / fig.get_figheight(), bottom=0.5 / fig.get_figheight())
    ymax = 1.05 * float(prof["p99_8"].max())
    use_um = prof["z_um"].notna().all()
    for i, s in enumerate(stacks):
        vmax = shared if scale == "shared" else vmaxes[s.series_name]
        idx = gallery_slice_indices(s.n_z)
        zum = s.z_um()
        im = None
        for j, (f, z) in enumerate(zip(GALLERY_FRACTIONS, idx, strict=True)):
            ax = fig.add_subplot(gs[i, j])
            im = ax.imshow(s.data[z], cmap="magma", vmin=0, vmax=vmax, interpolation="nearest")
            ax.set_xticks([]), ax.set_yticks([]), ax.grid(False)
            depth = f", {zum[z]:.0f} µm" if use_um else ""
            ax.set_title(f"{f:.0%}  (z={z}{depth})", fontsize=7.5)
            _scalebar(ax, s.pixel_size_um[1] if s.pixel_size_um else None, s.data.shape[1:])
            if j == 0:
                ax.set_ylabel(f"{s.concentration}\n{s.series_name}\n0–{vmax:.0f}\n{_short_unit(s.intensity_unit).lower()}",
                              rotation=0,
                              ha="right", va="center", fontsize=8.5, labelpad=8)
        cax = fig.add_subplot(gs[i, k])
        cb = fig.colorbar(im, cax=cax)  # unit and range are in the row label
        cb.ax.tick_params(labelsize=6.5)

        pax = fig.add_subplot(gs[i, k + 2])
        d = prof[prof["series_name"] == s.series_name].sort_values("z_index")
        x = d["z_um"] if use_um else d["z_index"]
        pax.plot(x, d["p99_8"], color=_color(s.concentration), lw=1.8)
        sel = d[d["z_index"].isin(idx)]
        pax.plot(sel["z_um"] if use_um else sel["z_index"], sel["p99_8"], "o", color=INK, markersize=5)
        pax.set_ylim(0, ymax)
        pax.set_ylabel("Slice p99.8 (counts)", fontsize=7.5)
        pax.tick_params(labelsize=7)
        if i == n - 1:
            pax.set_xlabel("Depth (µm, nominal)" if use_um else "z index", fontsize=8)
    mode = ("PER-IMAGE display (each row 0–p{p:g} of its own voxels; rows not comparable)" if scale != "shared"
            else "SHARED display 0–{v:.0f} for every row (largest per-stack p{p:g}; rows comparable)")
    fig.suptitle("Raw GFAP depth gallery (layout of Rose's first-staining figure): slices at 0/25/50/75/90% depth, "
                 "per-slice p99.8 at right (dots = shown slices)\n"
                 + mode.format(p=display_percentile, v=shared)
                 + ". Note: p99.8 near the detector ceiling flattens the surface peaks; see depth_profile_raw_mean.png",
                 fontsize=9.5)
    return _finish(fig, path, watermark)


def make_all(result: AnalysisResult, fig_dir: Path, display_percentile: float, zoom_max: float | None,
             watermark: str | None = None) -> list[Path]:
    out = [
        mip_grid(result, fig_dir / "raw_mip_grid_shared_scale.png", "shared", display_percentile, watermark),
        mip_grid(result, fig_dir / "raw_mip_grid_p998_display.png", "per_image", display_percentile, watermark),
        depth_profile(result, fig_dir / "depth_profile_raw_p998.png", "p99_8",
                      "Per-slice 99.8th percentile (raw counts)", watermark),
        depth_profile(result, fig_dir / "depth_profile_raw_mean.png", "mean", "Per-slice mean (raw counts)", watermark),
        depth_profile(result, fig_dir / "depth_profile_background_proxy.png", "bg_proxy_median",
                      "Background-proxy median (raw counts)", watermark,
                      note="Background PROXY per slice (not autofluorescence). Check background_summary.csv "
                           "bg_qc_flags before interpreting differences between concentrations."),
        background_comparison(result, fig_dir / "background_comparison.png", watermark),
        background_overlays(result, fig_dir / "background_ROI_overlays.png", display_percentile, watermark),
        threshold_sensitivity(result, fig_dir / "threshold_sensitivity.png", watermark),
        mip_with_histograms(result, fig_dir / "mip_with_histograms_rose_style.png", display_percentile, watermark),
        depth_slice_gallery(result, fig_dir / "depth_slice_gallery_per_image.png", "per_image", display_percentile,
                            watermark),
        depth_slice_gallery(result, fig_dir / "depth_slice_gallery_shared.png", "shared", display_percentile,
                            watermark),
    ]
    for src in ("mip_all_pixels", "stack_all_voxels", "mip_tissue", "stack_tissue"):
        out.append(histograms(result, fig_dir / f"histograms_{src}_shared_bins.png", src, zoom_max, watermark))
        out.append(survival(result, fig_dir / f"survival_{src}.png", src, watermark))
    return [p for p in out if p is not None]


# ---------------------------------------------------------------- review figures
def review_tissue_masks(stacks: list[GFAPStack], tissue: dict[str, np.ndarray], method: str, path: Path,
                        display_percentile: float = 99.8) -> Path:
    """MIPs (shared display) with the configured tissue mask outlined; outside is dimmed."""
    _style()
    mips = [max_projection(s.data) for s in stacks]
    vmax = shared_display_vmax(mips, display_percentile)
    fig, axes = plt.subplots(1, len(stacks), figsize=(3.3 * len(stacks), 3.8), squeeze=False, layout="constrained")
    for ax, s, mip in zip(axes[0], stacks, mips, strict=True):
        tm = tissue[s.series_name]
        ax.imshow(mip, cmap="magma", vmin=0, vmax=vmax, interpolation="nearest")
        if not tm.all():
            ax.imshow(np.where(tm, np.nan, 1.0), cmap="gray_r", alpha=0.55, vmin=0, vmax=1, interpolation="nearest")
            ax.contour(tm, levels=[0.5], colors=["#4fd1c5"], linewidths=1.0)
        ax.set_title(f"{s.concentration}  (raw MIP)\n{s.series_name}\ntissue = {tm.mean():.1%} of field",
                     fontsize=8.5)
        ax.set_xticks([]), ax.set_yticks([]), ax.grid(False)
        _scalebar(ax, s.pixel_size_um[1] if s.pixel_size_um else None, mip.shape)
    expl = ("whole field treated as tissue - confirm there is no empty space / mounting medium in any field"
            if method == "none" else "cyan outline = tissue; dimmed = excluded")
    fig.suptitle(f"REVIEW: tissue mask (method = {method}); {expl}\nShared display 0–{vmax:.0f} counts",
                 fontsize=9.5)
    return _finish(fig, path, "REVIEW - not a result")


def review_z_ranges(stacks: list[GFAPStack], tissue: dict[str, np.ndarray], configured: dict[str, tuple[int, int]],
                    auto: dict[str, tuple[int, int]], mode: str, path: Path, display_percentile: float = 99.8) -> Path:
    """Per stack: XZ max projection (max over Y) and per-slice profiles with the configured
    analysed range shaded and the auto-detected candidate marked."""
    _style()
    n = len(stacks)
    xz = [s.data.max(axis=1) for s in stacks]  # (Z, X)
    vmax = max(float(np.percentile(a, display_percentile)) for a in xz)
    fig, axes = plt.subplots(n, 2, figsize=(11, 2.6 * n), squeeze=False, layout="constrained",
                             gridspec_kw={"width_ratios": [1.3, 1]})
    for row, (s, img) in enumerate(zip(stacks, xz, strict=True)):
        ax_img, ax_p = axes[row]
        # true aspect would make a ~20-slice stack a hairline next to 2000+ px of X:
        # stretch z so the panel is ~1/3 as tall as wide, and say so in the title
        phys = (s.z_step_um / s.pixel_size_um[1]) if (s.z_step_um and s.pixel_size_um) else 1.0
        stretch = max(1.0, (img.shape[1] / 3.0) / (img.shape[0] * phys))
        ax_img.imshow(img, cmap="magma", vmin=0, vmax=vmax, aspect=phys * stretch, interpolation="nearest")
        z0, z1 = configured[s.series_name]
        a0, a1 = auto[s.series_name]
        for z in (z0 - 0.5, z1 + 0.5):
            ax_img.axhline(z, color="#4fd1c5", lw=1.2)
        for z in (a0 - 0.5, a1 + 0.5):
            ax_img.axhline(z, color="white", lw=0.9, ls=":")
        ax_img.set_ylabel("z index (0 = first acquired)")
        ax_img.set_xticks([]), ax_img.grid(False)
        ax_img.set_title(f"{s.concentration}: XZ max projection (max over Y), {s.series_name}; "
                         f"z stretched x{stretch:.1f}", loc="left", fontsize=8.5)

        tm = tissue[s.series_name]
        sl = [p[tm] if not tm.all() else p.ravel() for p in s.data]
        z = np.arange(s.n_z)
        ax_p.plot(z, [np.percentile(v, 99.8) for v in sl], color=INK, label="p99.8")
        ax_p.plot(z, [np.median(v) for v in sl], color=INK_2, ls="--", label="median")
        ax_p.axvspan(z0 - 0.5, z1 + 0.5, color="#4fd1c5", alpha=0.18, lw=0, label=f"configured ({mode})")
        ax_p.axvline(a0 - 0.5, color=FLAG_INK, lw=0.9, ls=":")
        ax_p.axvline(a1 + 0.5, color=FLAG_INK, lw=0.9, ls=":", label="auto candidate")
        ax_p.set_xlabel("z index")
        ax_p.set_ylabel("raw counts")
        ax_p.set_ylim(bottom=0)
        ax_p.set_title(f"configured {z0}–{z1}, auto candidate {a0}–{a1} of 0–{s.n_z - 1}", loc="left", fontsize=8.5)
        if row == 0:
            ax_p.legend(loc="upper left", bbox_to_anchor=(1.02, 1), fontsize=7.5)
    fig.suptitle("REVIEW: z-range. Cyan = analysed range as configured; dotted = auto-detected candidate "
                 "(unvalidated). Check that slices outside tissue are excluded and that nothing is cut from the tissue.",
                 fontsize=9.5)
    return _finish(fig, path, "REVIEW - not a result")
