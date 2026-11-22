# Validation status

_Last updated: 2026-11-20. No real microscopy data has been analysed yet._

## Verified now (synthetic data, `uv run pytest`, 90 tests)

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

## Confirmed from the real file (`gfap-conc inspect`, 2026-11-22)

`Sewon_GFAP_curve.lif` (plus a companion `Sewon_GFAP_curve.lifext`, not yet inspected):

- 8 images: `1-250`, `1-500`, `1-1000`, `1-2000` and `1-250_Lng` … `1-2000_Lng`. The
  suffix is written **`_Lng`** in the file (Rose wrote `_lng`). Name matching now ignores
  case, and the exact spelling is reported in `image_kind_evidence`.
- `_Lng` images differ in memory layout (CZYX vs ZCYX for raw) and reach 65535, with
  0.3–0.6% of AF647 voxels at 65535. This is consistent with deconvolution plus full-range
  rescaling.
- Channels (dye names from metadata): 0 = `DAPI`, 1 = `ALEXA 647`.
- Raw images: 2200 × 1500 px, z-step 1.0007 µm, 21/24/23/23 slices (about 20–23 µm, not
  50 µm), `Resolution` 16 bit.
- **Acquisition settings are identical for the four raw images** (from
  `acquisition_comparability_raw.csv`). Only z start/end, number of sections, stage XY
  position, z position and the setting name differ. Shared settings:
  - objective HC PL APO CS2 20x/0.75 **DRY**, zoom 2;
  - pixel 0.132 µm, field about 291 × 198 µm;
  - pinhole 1.0 AU;
  - scan speed 400, line accumulation 4, frame average 2;
  - 638 nm laser at 2.5% (AF647) and 405 nm at 2% (DAPI);
  - detectors HyD S 3 (AF647) and HyD S 2 (DAPI): SiPM in **PhotonCounting** mode, gain 2.5,
    offset 0.

  Absolute photon counts are therefore comparable at the acquisition level.
- **Four different stage positions** (X = 66.6 / 62.2 / 57.8 / 53.2 mm, about 4.4 mm
  apart). These are most likely four different sections, so region comparability (which
  part of the cord was imaged) must be confirmed with Rose.
- **Dry objective:** the axial step is nominal. Refractive-index mismatch between air and
  the mounting medium means the true depth into tissue is likely larger than the stage
  travel. Depth profiles are reported in nominal µm.
- **Nearly identical maxima** (AF647 1429/1421/1419/1418): only **one voxel** sits at each
  maximum and about 10⁻⁶ of voxels are within 2% of it, so there is **no hard clipping**.
  Soft compression of the brightest signal by photon-counting dead time cannot be excluded
  from these numbers. The upper tails of the survival curves will be checked, and the
  question raised with Rose.

### First review (2026-11-22, tissue mask `none`, z-range `full`)

- All four fields lie entirely inside tissue (no empty space, medium or section edge), so
  tissue mask `none` is appropriate. All fields share a layout: a dense bright band at the
  top and sparser processes below. Whether they are the same anatomical region must be
  confirmed with Rose.
- Per-slice p99.8 is high at **both** ends of every stack, with a dip in the middle (the
  pattern expected from penetration from both faces of a free-floating section).
  - First slices are at 67–89% of the peak, so the stacks start at or inside the tissue
    surface.
  - Last slices are at 32–73% of the peak with visible structure: approaching the lower
    surface, not outside it.
  - `full` is a defensible z-range.
- **The original `auto` rule failed on real data.** With `min + 0.5·(max − min)`, it
  proposed keeping only slices 0–6 (1:250) and 0–7 (1:2000), cutting at the mid-stack dip.
  It was not used. Fixed: the level is now measured from the background floor (lowest
  per-slice median), default 0.2, and only end slices can be trimmed. It remains
  unvalidated.
- The XZ review view is now stretched in z (factor printed) so that thin stacks are
  readable.

### Second review (2026-11-22, stretched XZ view, revised `auto` rule): z-range `full` confirmed

- No empty slices in any stack. End slices are darker in 1:500 (z22–23) and 1:1000 (z22),
  but their p99.8 is still ≥ ~400 counts against a background median of about 20, so they
  are tissue surfaces (partial volume), not medium. The revised `auto` rule independently
  keeps every slice (0–20, 0–23, 0–22, 0–22).
- **Comparability caveat:** the 1:250 stack spans 20.0 µm (stage Begin/End) against
  22.0–22.3 µm for the others. Its last slice is still at 73% of peak, so it probably
  stopped just short of the lower surface. Its last 20% window therefore holds the bright
  second surface, while the others include dim surface slices. `middle_over_last` is not
  strictly comparable across stacks; `middle_over_first` is the better retention measure
  (all stacks start at the upper surface, first slice at 67–89% of peak).

### First confirmed run (2026-11-22): what it changed in the pipeline

- Histograms used nominal 16-bit bins (256 counts wide) although all data are ≤ 1429
  counts, so the survival curves had only ~6 points. Bins now span the observed range by
  default (width 6 counts here).
- `LARGE_REGION` (> 85%) flagged all four stacks. GFAP fills only 2.3–4.4% of voxels at
  the shared threshold, so ~90% background is expected for this sparse marker. The limit is
  now 97%. No result value is affected. `RIGHT_TAIL` (tail ratio ≈ 5 in all stacks)
  remains, and is real: the proxy contains faint processes and out-of-focus light.
- `p99_8` lies within 5% of the near-identical maxima (~1355 vs ~1420) in every stack, so it
  sits in the range where any soft detector compression would act. p99.8-based depth
  ratios may understate the mid-depth dip; mean-based ratios are less exposed.
- **Re-run with 6-count bins (2026-11-22): soft detector compression is now strongly
  indicated.**
  - In all four stacks, the voxel-intensity density falls smoothly on a log scale up to
    ~1250 counts, then shows a **local bump at ~1300–1340 counts followed by a cliff to
    ~1420**. All four have the bump at the same place, although they are different sections.
  - A smooth biological distribution would not do this; compression of the brightest
    signal into a band below an effective ceiling would.
  - About 0.1% of stack voxels (above ~1250 counts) lie in this band. In the MIPs the share
    is much larger: MIP p99 is 1252–1327 and p99.8 is 1352–1364.
  - Consequence: maxima, MIP p99.8 and surface values of the per-slice p99.8 measure the
    ceiling, not the stain. Means, medians and the 100–1200-count part of the distribution
    are unaffected.
  - Confirmation requires a linearity test: the same field at two excitation powers or
    accumulation settings.

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
