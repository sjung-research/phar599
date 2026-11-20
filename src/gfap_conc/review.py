"""Mandatory human review of tissue masking and z-range before concentration comparison.

`gfap-conc review` produces QA figures and a checklist from the real images, but no
concentration-comparison numbers. The reviewer then records in the config:

    review:
      tissue_mask_confirmed: true
      z_range_confirmed: true
      reviewed_settings_sha256: <hash printed by `review`>

The hash covers the input file, channel, series and every tissue-mask / z-range setting.
If any of these changes after review, the confirmation no longer matches and `run`
refuses to compare concentrations (unless explicitly run as PRELIMINARY).
"""

from __future__ import annotations

import dataclasses
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from . import plots
from .config import ProjectConfig
from .model import GFAPStack, ImageKind, sort_by_dilution
from .pipeline import auto_z_range, compute_tissue_masks, resolve_z_ranges


def review_settings_sha256(cfg: ProjectConfig, lif_sha256: str) -> str:
    """Hash of everything the review approves: input, channel, series, mask and z-range settings."""
    payload = {
        "lif_sha256": lif_sha256,
        "gfap_channel": cfg.input.gfap_channel,
        "reference_channel": cfg.input.reference_channel,
        "bit_depth_override": cfg.input.bit_depth_override,
        "series": [[s.name, s.concentration, s.kind.value] for s in cfg.series],
        "tissue_mask": dataclasses.asdict(cfg.analysis.tissue_mask),
        "z_range": cfg.analysis.depth.z_range,
        "auto_metric": cfg.analysis.depth.auto_metric,
        "auto_rel_level": cfg.analysis.depth.auto_rel_level,
        "manual_z_ranges": {k: list(v) for k, v in sorted(cfg.analysis.depth.manual_z_ranges.items())},
    }
    return hashlib.sha256(yaml.safe_dump(payload, sort_keys=True).encode()).hexdigest()


def review_status(cfg: ProjectConfig, current_hash: str) -> tuple[bool, list[str]]:
    """(confirmed, reasons-if-not)."""
    r = cfg.review
    reasons = []
    if not r.tissue_mask_confirmed:
        reasons.append("review.tissue_mask_confirmed is not true")
    if not r.z_range_confirmed:
        reasons.append("review.z_range_confirmed is not true")
    if not r.reviewed_settings_sha256:
        reasons.append("review.reviewed_settings_sha256 is not set")
    elif r.reviewed_settings_sha256 != current_hash:
        reasons.append("input, series, tissue-mask or z-range settings changed since the review "
                       f"(config hash {r.reviewed_settings_sha256[:12]}…, current {current_hash[:12]}…)")
    return (not reasons), reasons


def build_review(
    stacks: list[GFAPStack],
    cfg: ProjectConfig,
    reference_images: dict[str, np.ndarray],
    out_dir: Path,
    settings_hash: str,
) -> list[Path]:
    raw = sort_by_dilution([s for s in stacks if s.kind is ImageKind.RAW])
    a = cfg.analysis
    tissue = compute_tissue_masks(raw, a, reference_images)
    configured, _ = resolve_z_ranges(raw, a, tissue)
    auto = {s.series_name: auto_z_range(s, tissue[s.series_name], a) for s in raw}

    rows = []
    for s in raw:
        tm = tissue[s.series_name]
        p998 = np.array([np.percentile(p[tm] if not tm.all() else p, 99.8) for p in s.data])
        peak = p998.max() if p998.max() > 0 else 1.0
        z0, z1 = configured[s.series_name]
        rows.append({
            "concentration": s.concentration,
            "series_name": s.series_name,
            "n_z": s.n_z,
            "z_step_um": s.z_step_um,
            "tissue_mask_method": a.tissue_mask.method,
            "tissue_fraction_of_field": float(tm.mean()),
            "z_range_mode": a.depth.z_range,
            "configured_first": z0,
            "configured_last": z1,
            "auto_candidate_first": auto[s.series_name][0],
            "auto_candidate_last": auto[s.series_name][1],
            "first_slice_p998_rel_to_peak": float(p998[0] / peak),
            "last_slice_p998_rel_to_peak": float(p998[-1] / peak),
            "peak_p998_slice": int(np.argmax(p998)),
        })
    table = pd.DataFrame(rows)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = [out_dir / "review_summary.csv"]
    table.to_csv(paths[0], index=False)
    fig_dir = out_dir / "figures"
    paths.append(plots.review_tissue_masks(raw, tissue, a.tissue_mask.method, fig_dir / "review_tissue_mask.png",
                                           a.display_percentile))
    paths.append(plots.review_z_ranges(raw, tissue, configured, auto, a.depth.z_range,
                                       fig_dir / "review_z_range.png", a.display_percentile))
    paths.append(_write_checklist(out_dir, table, cfg, settings_hash))
    return paths


def _write_checklist(out_dir: Path, table: pd.DataFrame, cfg: ProjectConfig, settings_hash: str) -> Path:
    hints = []
    for r in table.itertuples(index=False):
        if r.first_slice_p998_rel_to_peak < 0.5 or r.last_slice_p998_rel_to_peak < 0.5:
            hints.append(f"- `{r.series_name}`: an end slice has < 50% of the peak p99.8 "
                         f"(first {r.first_slice_p998_rel_to_peak:.2f}, last {r.last_slice_p998_rel_to_peak:.2f}); "
                         "check whether the stack extends outside the tissue.")
        if (r.auto_candidate_first, r.auto_candidate_last) != (r.configured_first, r.configured_last):
            hints.append(f"- `{r.series_name}`: auto candidate {r.auto_candidate_first}-{r.auto_candidate_last} "
                         f"differs from configured {r.configured_first}-{r.configured_last}.")
    text = f"""# Review checklist: tissue mask and z-range

Generated by `gfap-conc review`. This step must be completed on the **real** images
before `gfap-conc run` compares concentrations. No concentration-comparison numbers
were produced here.

Settings under review: tissue mask `{cfg.analysis.tissue_mask.method}` on `{cfg.analysis.tissue_mask.reference}`,
z-range `{cfg.analysis.depth.z_range}`.

## Look at

1. `figures/review_tissue_mask.png`
   - Does any field contain empty space, mounting medium, a tear or the section edge?
     If so, method `none` is wrong; choose `absolute` or `otsu` (optionally on the DAPI
     reference channel) and re-run `review`.
   - Does the mask follow the tissue without cutting into it?
   - Do the four fields come from comparable regions (grey/white matter, laminae, side)?
     Record this in `review.notes`. It cannot be fixed by masking.
2. `figures/review_z_range.png`
   - In the XZ view, where does tissue begin and end? Are slices outside the tissue
     included in the cyan range?
   - Is any tissue cut off by the configured range?
   - If trimming is needed, prefer `z_range: manual` with explicit `manual_z_ranges`
     chosen here, over the unvalidated `auto` rule.
   - Decide ranges from structure in the images, not from which concentration looks best.
3. `review_summary.csv`: numbers behind the figures.

## Automatic hints

{chr(10).join(hints) if hints else "_none_"}

## Record the decision in the config

When (and only when) both are satisfactory for the settings shown above:

```yaml
review:
  tissue_mask_confirmed: true
  z_range_confirmed: true
  reviewed_settings_sha256: {settings_hash}
  reviewed_by: <name>
  notes: "<what was checked; region comparability>"
```

Changing the input file, channel, series, tissue-mask or z-range settings changes this
hash, and the confirmation then no longer applies. Re-run `review` after any such change.
"""
    p = out_dir / "review_checklist.md"
    p.write_text(text, encoding="utf-8")
    return p
