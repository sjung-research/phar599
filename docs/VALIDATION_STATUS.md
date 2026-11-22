# Validation status

_Last updated: 2026-11-20. No real microscopy data has been analysed yet._

## Verified now (synthetic data, `uv run pytest`, 85 tests)

| Component | What the tests establish |
|---|---|
| Stack contract (`model.py`) | Non-3-D arrays (unresolved T/C) rejected; raw must be integer and within the declared bit depth; data are read-only; caller's array untouched |
| Raw vs `_lng` guard | Deconvolved stacks are rejected by intensity functions and absent from every output table |
| Statistics | Known means, medians, variance and percentiles; a single outlier moves `max` but not `p99_8`; no uint16 overflow |
| Analysis 1 | MIP taken over axis 0 (Z); a 2× brighter copy stays 2× in raw columns while the secondary normalised mean is identical; matched-depth windows match physical thickness |
| Analysis 2 | Known per-slice means and depths; within-slice variance kept separate from depth variation; non-overlapping windows; asymmetric profiles reported per side; z-range respected |
| Analysis 3 | XY-only exclusion radius; ROI coordinates map exactly to original pixels, out-of-bounds ROIs refused; pooled threshold equally weighted and deterministic; empty background region raises a QC warning; **shared-threshold bias for weak stains demonstrated** |
| Background QC | Tail asymmetry ≈ 1 for Poisson noise, > 2 when faint signal is mixed in; each flag rule triggers on its own; the synthetic weak stains are flagged and the "rises with dilution" warning fires; rows carry an explicit "not autofluorescence" note |
| Threshold sensitivity | Factor-1.0 rows equal the main background rows; thresholds are identical across concentrations at each factor; proxy fraction grows with T |
| Review gate | Hash changes with file, mask or z-range settings, but not with display settings; `run` refuses without confirmation, refuses after a reviewed setting changes, and marks `--allow-unreviewed` output PRELIMINARY; `review` writes no comparison tables; manual z-ranges validated |
| Origin separation | `data_origin` on every table; mixed origins refused; demo refuses non-`synthetic` paths; real runs refuse `synthetic` paths |
| Analysis 4 | Integer-aligned shared bins; `bincount` = `np.histogram`; totals equal pixel counts; partial last bin density-corrected; survival curves and quantiles exact |
| Pipeline | Identity columns (concentration, series name, channel) on every row; brightness order and mid-depth dip of synthetic data recovered; duplicate concentrations refused instead of pooled |
| LIF reader path | On small LIF files written to the documented v2 layout (`tests/lif_writer.py`): series listing, Z/C/Y/X axes, µm conversion, channel selection without mixing, bit depth, hardware-setting capture, comparability flags |
| CLI | `inspect` → `run` end to end on a test LIF; source file byte-identical afterwards; checksum mismatch aborts |

## Requires the real `.lif` before results can be trusted

Ordered roughly by when they will come up.

1. **Real LAS X file structure.** The LIF writer in the tests is ours, not LAS X output.
   We need to confirm that liffile opens the file, that series paths and names are as
   expected, and that the GFAP and DAPI channels sit at the indices we think. Check for
   unexpected axes: tile scans (`M`), time points, or a lambda dimension.
2. **Channel identity.** Confirm from `lif_xml_header.xml` (laser line about 633–653 nm,
   detection band, dye name) which `C` index is AF647. The code only refuses channels
   *named* DAPI/Hoechst and warns on a blue LUT. It cannot positively confirm AF647.
3. **Bit depth and intensity units.**
   - Does `Resolution` (e.g. 8/12/16) match the observed maximum value?
   - Was the detector in photon-counting or analog mode? This decides whether "counts"
     means photons.
   - Were accumulation or averaging used? Accumulation sums counts; averaging divides
     them.
4. **Hardware-setting attribute names.** `flatten_hardware` captures attributes
   generically from `ATLConfocalSettingDefinition`, `Detector`, `LaserLineSetting` and
   `MultiBand`. Real LAS X versions may name or nest these differently, and sequential
   scans have several setting blocks. Inspect `acquisition_comparability_raw.csv` by eye.
   Any difference in laser power, gain, zoom, pinhole or averaging between concentrations
   limits absolute comparisons.
5. **Voxel size convention.** Compare our `step = |Length| / (N − 1)` with the voxel size
   LAS X (or Fiji/Bio-Formats) displays for the same image.
6. **Tissue mask (enforced by `review`).** The default is `none` (whole field). In
   `review_tissue_mask.png`, check for empty space outside the section. If there is any,
   switch to `absolute` or `otsu`, possibly on the DAPI channel, re-run `review`, and
   confirm.
7. **Analysed z-range (enforced by `review`).** In `review_z_range.png`, check whether
   stacks start or end outside the tissue. Prefer explicit `manual_z_ranges` over the
   unvalidated `auto` rule. Decide from structure, before looking at concentration
   results.
8. **Foreground threshold and background proxy.**
   - Check that `pooled_otsu` gives a sensible GFAP mask on real images
     (`background_ROI_overlays.png`).
   - Read `threshold_sensitivity.png` and the background QC flags.
   - The QC limits (85% region, tail ratio 2, 10% threshold dependence) are first
     choices tuned on synthetic data. Check on the real images that they flag what an
     expert would flag, and no more.
   - If the proxy is flagged for some concentrations, prefer manual background ROIs chosen
     with Rose, or report background as not determinable.
9. **Saturation and zero clipping.** Review `qc_stack_summary.csv`. Saturation compresses
   high-concentration signal (one explanation for "same brightness"); zero clipping
   invalidates the background floor.
10. **Region comparability.** Are the four fields from equivalent spinal cord regions
    (grey vs white matter, laminae, ipsilateral vs contralateral)? This cannot be checked
    from pixel statistics. Ask Rose.
11. **Background-subtracted depth retention.** Raw middle/outer ratios are diluted by the
    background floor (see `docs/methods.md`). Once the background estimate is judged
    defensible on real data, decide whether to add `(metric − background)` ratios.
12. **Window sizes and stack thickness.** With the default 20% windows, k = round(0.2·n)
    slices per window: n = 30 gives 6, but n = 7 gives a single noisy slice per window.
    Fewer than 3 analysed slices produce a QC warning instead of a ratio. Confirm that the
    real slice counts give meaningful windows.

## Not implemented yet (deliberately deferred)

- `_lng` morphology QA (paired side-by-side crops, optional structural similarity). The
  `pair_of` field already exists in the config.
- Optional TIFF export of MIPs.
- Handling more than one field per concentration (replicates or ipsilateral/contralateral),
  which needs an explicit aggregation design.
- Morphology-aware histogram follow-up (ridge filters or local thickness) for the
  thin/thick-process hypothesis.
- Bio-Formats cross-reader check (byte-identical arrays from a second reader).
