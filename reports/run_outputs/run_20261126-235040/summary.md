# GFAP concentration curve - numeric summary (auto-generated)

Data origin: **lif**. Review status: **confirmed**.

Raw (non-deconvolved) series only. All intensities are unscaled source values (photon counts (HyD S PhotonCounting; line acc. 4, frame avg. 2)). No interpretation is contained in this file.

## Provenance

- `lif_path`: E:\2026NovemberGFAPlif\Sewon_GFAP_curve.lif
- `lif_sha256`: da134597fa839257395a3f5595f4c0a634910837337f518725bc278fc27041be
- `lif_size_bytes`: 2403467585
- `companion_files`: [{'path': 'E:\\2026NovemberGFAPlif\\Sewon_GFAP_curve.lifext', 'sha256': '534d60073c64663a78426e7efa98944217e695dc56be60deb896d9732db5e1a3', 'size_bytes': 384575069}]
- `config_path`: C:\Users\sjung3\Documents\phar599\config\analysis.yaml
- `config_sha256`: 8504d28050a2c1dcf8fd6b0416e37d5c69f9c61a8a2be22c41840b4c7c921a51
- foreground threshold: 341.0 (pooled_otsu (equal pixel sample per raw stack, shared by all concentrations))

## QC warnings

- Configured deconvolved series not loaded and excluded from all intensity analyses: 1-250_Lng, 1-500_Lng, 1-1000_Lng, 1-2000_Lng
- Raw stacks differ in number of z-slices: 1:250=21, 1:500=24, 1:1000=23, 1:2000=23
- Raw stacks differ in z-step (um): 1:250=1.000715, 1:500=1.0007152173913043, 1:1000=1.0007149999999998, 1:2000=1.0007149999999998
- 1-250 [1:250] background proxy RIGHT_TAIL: proxy intensities have a long upper tail, as expected if faint GFAP processes are included. Do not read this background value as autofluorescence.
- 1-500 [1:500] background proxy RIGHT_TAIL: proxy intensities have a long upper tail, as expected if faint GFAP processes are included. Do not read this background value as autofluorescence.
- 1-1000 [1:1000] background proxy RIGHT_TAIL: proxy intensities have a long upper tail, as expected if faint GFAP processes are included. Do not read this background value as autofluorescence.
- 1-2000 [1:2000] background proxy RIGHT_TAIL: proxy intensities have a long upper tail, as expected if faint GFAP processes are included. Do not read this background value as autofluorescence.
- Raw maxima are nearly identical across concentrations (1:250=1429, 1:500=1421, 1:1000=1419, 1:2000=1418) and far below the nominal ceiling (65535). Only a few voxels sit at each maximum, so there is no hard clipping; gradual (soft) compression of the brightest signal, e.g. photon-counting dead-time non-linearity, is NOT excluded by this check. saturation_fraction cannot detect either; compare the upper tails of the survival curves before comparing brightness.

## Stack QC

| concentration | series_name | shape_zyx | bit_depth | z_step_um | stack_max | stack_frac_saturated | stack_frac_zero | z_analysis_first | z_analysis_last |
|---|---|---|---|---|---|---|---|---|---|
| 1:250 | 1-250 | 21x1500x2200 | 16 | 1.001 | 1429 | 0 | 3.343e-05 | 0 | 20 |
| 1:500 | 1-500 | 24x1500x2200 | 16 | 1.001 | 1421 | 0 | 8.961e-05 | 0 | 23 |
| 1:1000 | 1-1000 | 23x1500x2200 | 16 | 1.001 | 1419 | 0 | 2.599e-05 | 0 | 22 |
| 1:2000 | 1-2000 | 23x1500x2200 | 16 | 1.001 | 1418 | 0 | 0.0001112 | 0 | 22 |

## 1. Max-intensity projection (raw)

`mean` is computed on unscaled raw MIP pixels. `display_vmax_within_image` (= per-image p99.8) is a display ceiling only; `normalized_mean_secondary` removes brightness differences by construction and must not be used to rank concentrations.

| concentration | projection | region | n_z_projected | mean | median | p99_8 | display_vmax_within_image | frac_saturated | normalized_mean_secondary |
|---|---|---|---|---|---|---|---|---|---|
| 1:250 | mip_full_stack | all_pixels | 21 | 219.5 | 110 | 1357 | 1357 | 0 | 0.1618 |
| 1:250 | mip_matched_depth_first | all_pixels | 21 | 219.5 | 110 | 1357 | 1357 | 0 | 0.1618 |
| 1:250 | mip_matched_depth_center | all_pixels | 21 | 219.5 | 110 | 1357 | 1357 | 0 | 0.1618 |
| 1:500 | mip_full_stack | all_pixels | 24 | 204.6 | 100 | 1352 | 1352 | 0 | 0.1513 |
| 1:500 | mip_matched_depth_first | all_pixels | 20 | 200.9 | 98 | 1351 | 1351 | 0 | 0.1487 |
| 1:500 | mip_matched_depth_center | all_pixels | 20 | 202.4 | 99 | 1349 | 1349 | 0 | 0.1501 |
| 1:1000 | mip_full_stack | all_pixels | 23 | 288.3 | 157 | 1364 | 1364 | 0 | 0.2114 |
| 1:1000 | mip_matched_depth_first | all_pixels | 21 | 287.7 | 156 | 1364 | 1364 | 0 | 0.2109 |
| 1:1000 | mip_matched_depth_center | all_pixels | 21 | 287.2 | 156 | 1364 | 1364 | 0 | 0.2105 |
| 1:2000 | mip_full_stack | all_pixels | 23 | 199.4 | 98 | 1354 | 1354 | 0 | 0.1473 |
| 1:2000 | mip_matched_depth_first | all_pixels | 21 | 197.5 | 96 | 1354 | 1354 | 0 | 0.1459 |
| 1:2000 | mip_matched_depth_center | all_pixels | 21 | 197.9 | 97 | 1353 | 1353 | 0 | 0.1463 |

## 2. Depth: middle vs outer windows

`first`/`last` = first/last acquired slices (orientation relative to the coverslip unknown).

| concentration | metric | window_n_slices | first_window_mean | middle_window_mean | last_window_mean | middle_over_first | middle_over_last | depth_cv |
|---|---|---|---|---|---|---|---|---|
| 1:250 | mean | 4 | 74.34 | 54.28 | 52.8 | 0.7301 | 1.028 | 0.2084 |
| 1:250 | median | 4 | 27.5 | 27.75 | 21.75 | 1.009 | 1.276 | 0.1575 |
| 1:250 | p99_8 | 4 | 1286 | 876.8 | 1008 | 0.6815 | 0.87 | 0.1662 |
| 1:500 | mean | 5 | 60.96 | 49.56 | 34.8 | 0.813 | 1.424 | 0.3152 |
| 1:500 | median | 5 | 22 | 25.8 | 14.6 | 1.173 | 1.767 | 0.2677 |
| 1:500 | p99_8 | 5 | 1197 | 812.8 | 802.2 | 0.6791 | 1.013 | 0.2464 |
| 1:1000 | mean | 5 | 86.68 | 75.53 | 60.1 | 0.8713 | 1.257 | 0.2452 |
| 1:1000 | median | 5 | 31.6 | 37.2 | 23.2 | 1.177 | 1.603 | 0.2427 |
| 1:1000 | p99_8 | 5 | 1282 | 1108 | 1079 | 0.8643 | 1.027 | 0.116 |
| 1:2000 | mean | 5 | 64.52 | 49.24 | 43.39 | 0.7632 | 1.135 | 0.2383 |
| 1:2000 | median | 5 | 23.6 | 24.8 | 16.8 | 1.051 | 1.476 | 0.2017 |
| 1:2000 | p99_8 | 5 | 1265 | 850.6 | 968.6 | 0.6724 | 0.8782 | 0.1753 |

## 3. Background (proxy unless stated)

A background **proxy** is low-intensity tissue away from the shared-threshold foreground. It is not a measurement of autofluorescence. A proxy that is higher for one concentration may reflect faint GFAP below the threshold, not tissue autofluorescence. Check `bg_qc_flags` and the threshold-sensitivity table first. `bg_robust_sd` is within-image spread (1.4826 × MAD), not replicate uncertainty.

| concentration | background_method | bg_median | bg_robust_sd | bg_fraction_of_reference | bg_tail_asymmetry | bg_threshold_dependence | bg_qc_flags | fg_median | signal_to_background | contrast_to_noise |
|---|---|---|---|---|---|---|---|---|---|---|
| 1:250 | proxy_stack_3d | 24 | 19.27 | 0.9174 | 5.294 | 0.08333 | RIGHT_TAIL | 542 | 22.58 | 26.88 |
| 1:250 | proxy_mip_2d (inflated by max over z) | 59 | 41.51 | 0.5618 | 3.098 | – | visual QA only | 603 | 10.22 | 13.1 |
| 1:500 | proxy_stack_3d | 21 | 16.31 | 0.9346 | 5.25 | 0.04762 | RIGHT_TAIL | 541 | 25.76 | 31.89 |
| 1:500 | proxy_mip_2d (inflated by max over z) | 58 | 40.03 | 0.6 | 3.125 | – | visual QA only | 606 | 10.45 | 13.69 |
| 1:1000 | proxy_stack_3d | 28 | 23.72 | 0.8885 | 4.952 | 0.07143 | RIGHT_TAIL | 558 | 19.93 | 22.34 |
| 1:1000 | proxy_mip_2d (inflated by max over z) | 70 | 45.96 | 0.4449 | 2.667 | – | visual QA only | 640 | 9.143 | 12.4 |
| 1:2000 | proxy_stack_3d | 21 | 17.79 | 0.9314 | 5.438 | 0.04762 | RIGHT_TAIL | 542 | 25.81 | 29.28 |
| 1:2000 | proxy_mip_2d (inflated by max over z) | 58 | 41.51 | 0.6167 | 3.049 | – | visual QA only | 599 | 10.33 | 13.03 |

### Threshold sensitivity (diagnostic)

Same exclusion rule re-run at multiples of the one shared threshold, identically for all concentrations. The analysis itself uses factor 1.0.

| concentration | threshold_factor | foreground_threshold | fg_fraction_of_tissue | bg_fraction_of_reference | bg_median | bg_tail_asymmetry | signal_to_background |
|---|---|---|---|---|---|---|---|
| 1:250 | 0.5 | 170.5 | 0.06953 | 0.8233 | 22 | 3.375 | 13.36 |
| 1:250 | 0.75 | 255.8 | 0.04173 | 0.8852 | 23 | 4.625 | 18.43 |
| 1:250 | 1 | 341 | 0.02846 | 0.9174 | 24 | 5.294 | 22.58 |
| 1:250 | 1.25 | 426.2 | 0.02067 | 0.9371 | 25 | 5.667 | 26 |
| 1:250 | 1.5 | 511.5 | 0.01562 | 0.9507 | 25 | 6.333 | 29.96 |
| 1:250 | 2 | 682 | 0.009461 | 0.9682 | 26 | 7 | 35.23 |
| 1:500 | 0.5 | 170.5 | 0.0569 | 0.8559 | 19 | 3.786 | 15.42 |
| 1:500 | 0.75 | 255.8 | 0.03411 | 0.908 | 20 | 4.733 | 21.15 |
| 1:500 | 1 | 341 | 0.02325 | 0.9346 | 21 | 5.25 | 25.76 |
| 1:500 | 1.25 | 426.2 | 0.01686 | 0.9508 | 21 | 6 | 30.86 |
| 1:500 | 1.5 | 511.5 | 0.01273 | 0.9616 | 21 | 6.625 | 35.43 |
| 1:500 | 2 | 682 | 0.00765 | 0.9753 | 22 | 7.118 | 41.27 |
| 1:1000 | 0.5 | 170.5 | 0.1017 | 0.7737 | 24 | 3.278 | 12.75 |
| 1:1000 | 0.75 | 255.8 | 0.06342 | 0.8481 | 27 | 4.15 | 16.22 |
| 1:1000 | 1 | 341 | 0.04434 | 0.8885 | 28 | 4.952 | 19.93 |
| 1:1000 | 1.25 | 426.2 | 0.03285 | 0.9142 | 29 | 5.545 | 23 |
| 1:1000 | 1.5 | 511.5 | 0.02532 | 0.9318 | 30 | 5.957 | 25.47 |
| 1:1000 | 2 | 682 | 0.01582 | 0.9549 | 31 | 6.875 | 30 |
| 1:2000 | 0.5 | 170.5 | 0.06027 | 0.8513 | 19 | 3.857 | 15.37 |
| 1:2000 | 0.75 | 255.8 | 0.03601 | 0.9044 | 20 | 4.867 | 21.1 |
| 1:2000 | 1 | 341 | 0.02443 | 0.9314 | 21 | 5.438 | 25.81 |
| 1:2000 | 1.25 | 426.2 | 0.0177 | 0.9484 | 21 | 6.188 | 31 |
| 1:2000 | 1.5 | 511.5 | 0.01341 | 0.9596 | 22 | 6.412 | 34 |
| 1:2000 | 2 | 682 | 0.008116 | 0.9738 | 22 | 7.353 | 41.73 |

## 4. Histogram quantiles (exact, from integer level counts)

| concentration | source | n_pixels | q0.5 | q0.9 | q0.99 | q0.998 |
|---|---|---|---|---|---|---|
| 1:250 | mip_all_pixels | 3300000 | 110 | 584 | 1278 | 1357 |
| 1:250 | stack_all_voxels | 69300000 | 27 | 126 | 662 | 1146 |
| 1:500 | mip_all_pixels | 3300000 | 100 | 543 | 1252 | 1352 |
| 1:500 | stack_all_voxels | 79200000 | 22 | 106 | 590 | 1082 |
| 1:1000 | mip_all_pixels | 3300000 | 157 | 785 | 1327 | 1364 |
| 1:1000 | stack_all_voxels | 75900000 | 33 | 173 | 849 | 1243 |
| 1:2000 | mip_all_pixels | 3300000 | 98 | 521 | 1254 | 1354 |
| 1:2000 | stack_all_voxels | 75900000 | 23 | 112 | 609 | 1113 |

## Figures

- `figures/raw_mip_grid_shared_scale.png`
- `figures/raw_mip_grid_p998_display.png`
- `figures/depth_profile_raw_p998.png`
- `figures/depth_profile_raw_mean.png`
- `figures/depth_profile_background_proxy.png`
- `figures/background_comparison.png`
- `figures/background_ROI_overlays.png`
- `figures/threshold_sensitivity.png`
- `figures/mip_with_histograms_rose_style.png`
- `figures/depth_slice_gallery_per_image.png`
- `figures/depth_slice_gallery_shared.png`
- `figures/histograms_mip_all_pixels_shared_bins.png`
- `figures/survival_mip_all_pixels.png`
- `figures/histograms_stack_all_voxels_shared_bins.png`
- `figures/survival_stack_all_voxels.png`

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
