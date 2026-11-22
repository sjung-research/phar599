"""Run manifest (provenance) and an auto-generated numeric summary.

``summary.md`` contains numbers and QC flags only. Interpretation is a separate,
human-written step; the file ends with the template from the analysis spec.
"""

from __future__ import annotations

import json
import platform
import subprocess
import sys
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path

import numpy as np
import pandas as pd

from .pipeline import AnalysisResult

PACKAGES = ("gfap-conc", "numpy", "scipy", "pandas", "matplotlib", "scikit-image", "tifffile", "liffile", "pyyaml")


def _git(*args: str) -> str | None:
    try:
        return subprocess.run(["git", *args], capture_output=True, text=True, check=True, encoding="utf-8", errors="replace",
                              cwd=Path(__file__).resolve().parent).stdout.strip()
    except Exception:
        return None


def write_manifest(out_dir: Path, *, config: dict, config_sha256: str, provenance: dict,
                   result: AnalysisResult, synthetic: bool, review: dict) -> Path:
    versions = {}
    for p in PACKAGES:
        try:
            versions[p] = metadata.version(p)
        except metadata.PackageNotFoundError:
            versions[p] = None
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "synthetic_data": synthetic,
        "data_origin": result.params.get("data_origin"),
        "review_status": review,
        "git_commit": _git("rev-parse", "HEAD"),
        "git_dirty": bool(_git("status", "--porcelain")),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "packages": versions,
        "input": provenance,
        "config_sha256": config_sha256,
        "config_resolved": config,
        "analysis_params": result.params,
        "raw_series": [s.identity() | {"source": s.source} for s in result.raw_stacks],
        "excluded_series": [s.identity() | {"source": s.source} for s in result.excluded_stacks],
        "warnings": result.warnings,
    }
    p = out_dir / "run_manifest.json"
    p.write_text(json.dumps(manifest, indent=2, default=str), encoding="utf-8")
    return p


def md_table(df: pd.DataFrame, columns: list[str]) -> str:
    cols = [c for c in columns if c in df.columns]
    if df.empty or not cols:
        return "_(no rows)_\n"

    def fmt(v) -> str:
        if isinstance(v, (float, np.floating)):
            if np.isnan(v):
                return "–"
            return f"{v:.4g}"
        return str(v)

    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    lines += ["| " + " | ".join(fmt(v) for v in row) + " |" for row in df[cols].itertuples(index=False)]
    return "\n".join(lines) + "\n"


INTERPRETATION_TEMPLATE = """\
## Interpretation (to be written by a person after reviewing QC)

For each concentration, address:

- Raw MIP mean and p99.8 (units; display assumptions)
- Middle vs outer-depth signal retention
- Background / signal-to-background (state whether proxy or validated ROI)
- Histogram shape and *plausible* morphological interpretation (histograms alone
  cannot establish process thickness)
- Acquisition / dtype / ROI comparability and sample count (n = 1 field per
  concentration: no inferential statistics; pixels and slices are not replicates)
- Conclusion and confidence: strongest evidence, trade-offs, unknowns.
  "No defensible optimal concentration from these images" is an acceptable outcome.
"""


def write_summary(out_dir: Path, result: AnalysisResult, *, provenance: dict, config_sha256: str,
                  synthetic: bool, figures: list[Path], review: dict) -> Path:
    t = result.tables
    parts = ["# GFAP concentration curve - numeric summary (auto-generated)\n"]
    parts.append(f"Data origin: **{result.params.get('data_origin')}**. Review status: **{review['status']}**.\n")
    if synthetic:
        parts.append("> **SYNTHETIC DATA.** Generated to exercise the pipeline. These numbers are not experimental "
                     "results and must not be reported.\n")
    elif review["status"] != "confirmed":
        parts.append("> **PRELIMINARY.** Tissue mask and/or z-range were not confirmed by review for these settings "
                     "(" + "; ".join(review.get("reasons", [])) + "). Do not compare concentrations from this run.\n")
    parts.append("Raw (non-deconvolved) series only. All intensities are unscaled source values "
                 f"({result.raw_stacks[0].intensity_unit}). No interpretation is contained in this file.\n")
    parts.append("## Provenance\n")
    for k, v in provenance.items():
        parts.append(f"- `{k}`: {v}")
    parts.append(f"- `config_sha256`: {config_sha256}")
    parts.append(f"- foreground threshold: {result.params.get('foreground_threshold')} "
                 f"({result.params.get('foreground_threshold_source')})\n")

    parts.append("## QC warnings\n")
    parts.append("\n".join(f"- {w}" for w in result.warnings) + "\n" if result.warnings else "_none_\n")
    parts.append("## Stack QC\n")
    parts.append(md_table(t["qc_stack_summary"], ["concentration", "series_name", "shape_zyx", "bit_depth",
                                                   "z_step_um", "stack_max", "stack_frac_saturated",
                                                   "stack_frac_zero", "z_analysis_first", "z_analysis_last"]))

    parts.append("## 1. Max-intensity projection (raw)\n")
    parts.append("`mean` is computed on unscaled raw MIP pixels. `display_vmax_within_image` (= per-image "
                 "p99.8) is a display ceiling only; `normalized_mean_secondary` removes brightness differences "
                 "by construction and must not be used to rank concentrations.\n")
    mp = t["max_projection_summary"]
    parts.append(md_table(mp, ["concentration", "projection", "region", "n_z_projected", "mean", "median",
                               "p99_8", "display_vmax_within_image", "frac_saturated",
                               "normalized_mean_secondary"]))

    parts.append("## 2. Depth: middle vs outer windows\n")
    parts.append("`first`/`last` = first/last acquired slices (orientation relative to the coverslip unknown).\n")
    parts.append(md_table(t["middle_vs_outer_summary"], ["concentration", "metric", "window_n_slices",
                                                          "first_window_mean", "middle_window_mean",
                                                          "last_window_mean", "middle_over_first",
                                                          "middle_over_last", "depth_cv"]))

    parts.append("## 3. Background (proxy unless stated)\n")
    parts.append("A background **proxy** is low-intensity tissue away from the shared-threshold foreground. "
                 "It is not a measurement of autofluorescence. A proxy that is higher for one concentration may "
                 "reflect faint GFAP below the threshold, not tissue autofluorescence. Check `bg_qc_flags` "
                 "and the threshold-sensitivity table first. `bg_robust_sd` is within-image spread "
                 "(1.4826 × MAD), not replicate uncertainty.\n")
    parts.append(md_table(t["background_summary"], ["concentration", "background_method", "bg_median",
                                                     "bg_robust_sd", "bg_fraction_of_reference", "bg_tail_asymmetry",
                                                     "bg_threshold_dependence", "bg_qc_flags", "fg_median",
                                                     "signal_to_background", "contrast_to_noise"]))
    sens = t.get("threshold_sensitivity")
    if sens is not None and not sens.empty:
        parts.append("### Threshold sensitivity (diagnostic)\n")
        parts.append("Same exclusion rule re-run at multiples of the one shared threshold, identically for all "
                     "concentrations. The analysis itself uses factor 1.0.\n")
        parts.append(md_table(sens, ["concentration", "threshold_factor", "foreground_threshold",
                                     "fg_fraction_of_tissue", "bg_fraction_of_reference", "bg_median",
                                     "bg_tail_asymmetry", "signal_to_background"]))

    parts.append("## 4. Histogram quantiles (exact, from integer level counts)\n")
    parts.append(md_table(t["histogram_quantiles"], ["concentration", "source", "n_pixels", "q0.5", "q0.9",
                                                      "q0.99", "q0.998"]))

    parts.append("## Figures\n")
    parts.append("\n".join(f"- `figures/{p.name}`" for p in figures) + "\n")
    parts.append(INTERPRETATION_TEMPLATE)
    p = out_dir / "summary.md"
    p.write_text("\n".join(parts), encoding="utf-8")
    return p
