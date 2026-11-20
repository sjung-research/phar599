# Methods: calculation definitions and assumptions

All quantities are computed on **raw, unscaled pixel values** of the GFAP (AF647)
channel of the **raw** image series, in source units ("raw counts"). Whether a count is
a photon count (HyD counting mode) or an analog digital level depends on the
acquisition mode, which must be read from metadata. Statistics are accumulated in float64.
Percentiles use numpy `method="linear"` (type 7), so `p99_8` may be non-integer.

Notation: `I[z, y, x]` is the raw stack, z = 0 is the **first acquired slice**, Z is the
number of slices, and `sat = 2**bits − 1`, where `bits` is the true acquisition bit depth
from the channel's `Resolution` metadata.

## 0. Data handling

- **Axes:** liffile labels are used to select one channel along `C`. A singleton `T` is
  required. Mosaic, RGB, spectral and time axes raise an error. The result is transposed
  to explicit `(Z, Y, X)`. There are no reshapes by assumed order.
- **Physical sizes:** `step = |Length| / (N − 1)`, converted from the dimension's `Unit`
  (LAS X stores metres). The sign of the Z `Length` is recorded as `z_direction`.
- **Immutability:** stacks are read-only views. The `.lif` is opened with mode `r` and
  verified by SHA-256 before every run.
- **Exclusion:** stacks of kind `deconvolved` (`_lng`) raise `DeconvolvedIntensityError` in
  every intensity function, and the pipeline lists them in `qc_warnings.txt`.
- **One stack per concentration:** more than one raw stack for a concentration is an
  error. Replicate or field aggregation must be designed explicitly, not pooled silently.

## 1. Max-intensity projection

- `MIP[y, x] = max_z I[z, y, x]` (full stack; `projection = mip_full_stack`).
- `mean, median, std, variance, p95, p99_8, min, max` are computed over MIP pixels
  (`region = all_pixels`, plus `tissue_mask` if a partial tissue mask is used).
- `frac_saturated = #(MIP ≥ sat) / n`, and `frac_zero = #(MIP == 0) / n`.
- `display_vmax_within_image = percentile(MIP, 99.8)`: Rose's "normalization max", used
  **only** for display.
- `normalized_mean_secondary = mean / display_vmax_within_image`. **Secondary:** it equals
  the same value for an image and its 2× brighter copy (see `tests/test_projection.py`).
- **Matched-depth MIPs** (`mip_matched_depth_first` / `_center`) project over the same
  physical thickness in every stack, namely that of the thinnest analysed stack. They are
  emitted only when stacks differ in depth, because a MIP's mean grows with the number of
  slices projected.
- **Figures:** the shared scale uses `vmax = max over concentrations of the per-image
  p99.8`, with `vmin = 0`. The per-image figure uses each image's own p99.8, and its title
  warns that brightness is not comparable there.

## 2. Depth profile

- Per slice z: the same statistics as above, over the slice (or tissue mask), plus:
  - `z_um = z · z_step`;
  - `depth_fraction_stack = z / (Z − 1)`;
  - `depth_fraction_analysis`, the same within the analysed z-range.
- `variance` is the **within-slice spatial** pixel variance. Variation **across depth** is
  summarised separately as `depth_cv = SD(profile) / mean(profile)` and
  `min_over_max`. These are different quantities.
- `foreground_fraction` is the fraction of slice pixels ≥ the shared foreground threshold.
- `bg_proxy_median` is the median of background-proxy pixels in that slice (section 3).
- **Windows:** for n analysed slices, `k = max(1, round(0.2 n))`:
  - first = `[0, k)`;
  - last = `[n − k, n)`;
  - middle = k slices starting at `(n − k) // 2`.

  At least `3k ≤ n` slices are required. Each window value is the mean of the per-slice
  metric (`mean`, `median`, `p99_8`). The output reports `middle_over_first`,
  `middle_over_last` and `middle_over_outer_mean`. First and last are reported separately,
  because orientation relative to the coverslip is unknown and profiles are rarely
  symmetric.
- **Analysed z-range:** `full` by default. `auto` keeps slices between the first and last
  crossing of `min + 0.5·(max − min)` of the per-slice `p99_8`. It is meant to trim slices
  outside the tissue and is **unvalidated**.
- Nominal z is not corrected for refractive-index mismatch.
- **Pitfall: ratios on raw values include the background floor.** An additive background
  that does not attenuate with depth pulls `middle_over_*` toward 1, so a **dimmer stain can
  look more uniform with depth** than a bright one even when the GFAP-specific signal
  attenuates identically (visible in the synthetic demo). Compare `p99_8`-based ratios,
  which are less affected, and treat retention differences between concentrations with
  this in mind. A background-subtracted retention variant is deferred until a defensible
  background estimate exists for the real data.

## 3. Background

- **Shared foreground threshold** `T` (raw counts) is either fixed in the config or
  computed as `pooled_otsu`: Otsu's threshold on an equal-sized random sample of 200,000
  tissue voxels from each raw stack (seed 0). It is one value, used for every concentration.
- **Background proxy (primary, `proxy_stack_3d`):** tissue voxels with `I < T` that are
  more than `exclusion_radius_px` (default 5) from any foreground voxel in the same slice.
  Dilation is in XY only. The summary reports `bg_median`, `bg_mean`,
  `bg_robust_sd = 1.4826·MAD`, `bg_p05`, `bg_p25`, `bg_frac_zero` and the fraction of tissue
  voxels used.
- **`proxy_mip_2d`:** the same rule applied to the MIP. The maximum over z inflates the
  background floor, so this is labelled as such and kept for visual QA.
- **`manual_roi_stack`:** half-open rectangles `[y0, y1) × [x0, x1)` in original pixel
  coordinates, applied to every slice. Out-of-bounds ROIs are an error; nothing is clipped
  or resized.
- **Derived (proxy-based) metrics:**
  - `signal_to_background = fg_median / bg_median`;
  - `contrast_to_noise = (fg_median − bg_median) / bg_robust_sd`;

  where `fg` = tissue voxels ≥ T.
- **Known bias:** with one shared threshold, faint real processes in a weak stain fall
  below T and enter the background proxy. This raises its median and lowers the
  signal-to-background ratio for dim concentrations, as shown by
  `test_background_rows_labelled_as_proxy`. Check `figures/background_ROI_overlays.png`,
  compare `bg_p05`/`bg_p25` with `bg_median`, and treat signal-to-background differences
  between concentrations cautiously.
- A low-intensity proxy is **not** a measurement of autofluorescence. Perfusion
  autofluorescence, non-specific binding and detector offset cannot be distinguished
  without no-primary controls.

## 4. Histograms

- **Bins:** `width = ceil((sat + 1) / n_bins)`, with integer edges `0, width, 2·width, …`
  until `sat` is covered. Every bin spans the same number of integer levels; the last may
  hold fewer, recorded in `n_levels_in_bin`. The same edges are used for all
  concentrations.
- **Counting:** `np.bincount` on exact integer values. Totals are asserted to equal the
  pixel count.
- **Columns:**
  - `fraction = count / n`;
  - `density_per_level = fraction / n_levels_in_bin`, which is what is plotted (log y-axis).
- **Sources, kept separate because they answer different questions:**
  - `mip_all_pixels`: projected structure;
  - `stack_all_voxels`: every voxel in the analysed z-range;
  - `*_tissue` variants when a tissue mask is used.
- **Survival:** `fraction_at_or_above(T)` for each bin edge T. This is the fraction of
  pixels above any threshold, so it avoids committing to a single threshold.
- **Quantiles:** `q0.5 … q0.998` are lower inverted-CDF quantiles from the level counts.
  They can differ slightly from the interpolated `p99_8` in other tables.
- Rose's thin-vs-thick-process hypothesis is a morphological interpretation. Histogram
  shape alone cannot establish process thickness.

## Statistical stance

There is one field per concentration and no biological replicates. Results are descriptive.
No p-values are computed, and pixels or slices are never treated as independent samples.
Differences between concentrations are confounded with field and section position
(e.g. white vs grey matter), so region comparability must be checked before
interpretation.
