# Progress summary: Mouse anti-GFAP concentration curve (input for writing a progress report)

**Student:** Sewon (PHAR 599 D1/D2, McGill University)
**Lab:** Prof. Alfredo Ribeiro-da-Silva
**Direct mentor:** Rose Rodrigues
**Period covered:** late October to 27 November 2026, with the analysis work done 18–23 November 2026.

---

## 0. Instructions for the AI writing the progress report

- **Use only the facts in this file.** Do not invent numbers, methods, dates or conclusions.
  If something the report needs is missing, write "[to be added]" rather than guessing.
- **Keep the scientific caution.** There is **one image (field) per concentration**, from
  four different tissue sections, so:
  - no p-values or claims of statistical significance;
  - no claim of a definitive optimal concentration;
  - findings are "descriptive" and "preliminary".
- **Separate what the student did from what others did.**
  - Rose prepared the serial dilution and acquired the images.
  - The student did the Day-2 staining (secondary, DAPI, washes) and the quantitative
    analysis.
- **Disclose AI assistance.** The analysis code was written with an AI coding assistant
  (Claude) under the student's direction. The student defined the requirements, ran
  everything on the real data, and made the review decisions (tissue mask, z-range).
  The report should say so.
- **Language:** write the report in English unless asked otherwise.
- **Figures:** the figures listed in section 7 exist as PNG files. If they are not attached,
  describe them in words as section 7 does.

---

## 1. Background and aim

- **Long-term lab goal:** build a reproducible imaging-to-analysis pipeline to visualise,
  annotate and eventually segment individual protoplasmic astrocytes in mouse spinal cord.
  This would use GFAP immunostaining, sparse tdTomato labelling and Sox9, with napari-based
  manual correction of threshold pre-masks. Good, consistent GFAP staining is a prerequisite.
- **Earlier work (context; from lab communications):**
  - Mouse, rabbit and chicken anti-GFAP primaries were compared, all with Alexa Fluor 647
    secondaries.
  - All gave signal, but the middle of the tissue stacks was dimmer than the surfaces
    ("penetration issue"); rabbit was worst.
  - Rose analysed this with per-slice 99.8th-percentile depth profiles.
  - A new Mouse anti-GFAP lot matching Rose's original successful antibody was ordered.
- **Current experiment (Aim 1: staining optimisation):** a 2-fold serial dilution of the new
  Mouse anti-GFAP primary at **1:250, 1:500, 1:1000 and 1:2000**, with Donkey anti-Mouse
  AF647 secondary and DAPI, on fixed free-floating mouse spinal cord sections (nominally
  50 µm).
  - Rose made the dilutions.
  - The student completed Day 2 (secondary, DAPI, washes).
  - Rose imaged the sections.
- **Why the analysis was needed:** Rose observed that this aliquot "looks exactly the same in
  all concentrations", unlike a colleague's (Emma's) aliquot, which changed with
  concentration. So no final concentration could be chosen.
- **Rose's four requested analyses:**
  1. average bit value of max-projected images, using a 99.8th-percentile "normalization
     max";
  2. intensity change through the z-stack ("how much darker is it in the middle?");
  3. differences in background fluorescence;
  4. intensity histogram distribution (hypothesis: thin processes give a low-intensity
     peak, thick processes skew high).
- **Rose's constraint:** the deconvolved images (`_Lng` suffix) are rescaled to the full
  16-bit range and must not be used to decide which concentration is best.

## 2. What was done

1. **Planning (Nov 20).** Read the project documents and designed a workflow:
   - inputs, calculations, QC checks, outputs;
   - pitfalls: per-image normalisation erasing real differences, detector saturation,
     background estimation, depth alignment, field-to-field variability.
2. **Software.** Built a Python package, `gfap-conc`, with:
   - modular code for each of Rose's four analyses;
   - reproducible environment management (uv, locked dependencies);
   - 90 automated tests on synthetic data;
   - documentation of every calculation.

   It is in a private GitHub repository (`sjung-research/phar599`).
3. **Built-in safeguards:**
   - the raw `.lif` is read-only and verified by checksum on every run;
   - deconvolved `_Lng` images are blocked from all intensity calculations;
   - intensities are never rescaled before comparison (the 99.8th percentile is kept for
     display only);
   - one shared foreground threshold is used for all concentrations;
   - background is labelled a "proxy";
   - a mandatory human review of the tissue mask and z-range must be passed before any
     comparison is run;
   - every run writes a provenance record.
4. **Real-data inspection (Nov 22).** Read the Leica file's metadata, which confirmed:
   - series names, the AF647 channel, voxel sizes and detector mode;
   - that **all intensity-relevant acquisition settings were identical** across the four
     images.
5. **Review (Nov 22).**
   - The student checked tissue-mask and z-range QA figures and confirmed that all four
     fields lie fully within tissue and that no slices are empty.
   - Problems found and fixed: an automatic z-range rule failed on the real data, and the
     side-view figure was unreadable.
6. **Analysis runs (Nov 22).** Ran all four analyses. Two issues revealed by the real data
   were fixed: histogram bins were too coarse, and one QC flag was miscalibrated. Re-run.
7. **Reproducibility (Nov 23).**
   - The same analysis ran on a second computer (the lab desktop) with a checksum-verified
     copy of the data, giving the same outputs and warnings.
   - A byte-level comparison of the output tables between the two computers was suggested
     but not yet reported.
8. **Documentation (Nov 23).**
   - An English results report (with 7 figures, also as .docx) for Rose;
   - methods and validation documents;
   - a session log;
   - a snapshot of the full run output, stored in the repository.

## 3. Data and acquisition (from file metadata)

| Item | Value |
|---|---|
| File | `Sewon_GFAP_curve.lif` (2.40 GB) + `.lifext` companion (0.38 GB, not needed) |
| Images | 4 raw z-stacks (`1-250`, `1-500`, `1-1000`, `1-2000`) + 4 deconvolved (`_Lng`, excluded) |
| Channels | 0 = DAPI (405 nm laser, 2%); **1 = Alexa 647 / GFAP** (638 nm laser, 2.5%) |
| Microscope settings (identical for all four) | 20x/0.75 dry objective, zoom 2; HyD S detectors in **photon-counting** mode, gain 2.5; line accumulation 4, frame average 2; pinhole 1 AU |
| Image size | 2200 × 1500 px at 0.132 µm (about 291 × 198 µm field) |
| z | 1.0 µm steps; 21, 24, 23 and 23 slices (about 20–22 µm, nominal) |
| Sections | Four different stage positions, i.e. four different sections (one field each) |

## 4. Results

All values are raw photon counts from the confirmed analysis run.

### 4.1 Max-projection intensity (Rose's "average bit value")

| Dilution | MIP mean | MIP median | MIP 99.8th percentile | GFAP-positive voxels (≥ shared threshold) |
|---|---|---|---|---|
| 1:250 | 219.5 | 110 | 1357 | 2.8% |
| 1:500 | 204.6 | 100 | 1352 | 2.3% |
| 1:1000 | 288.3 | 157 | 1364 | 4.4% |
| 1:2000 | 199.4 | 98 | 1354 | 2.4% |

- **No dose dependence.** 8× more primary (1:2000 → 1:250) raised MIP mean by only 10%.
- 1:1000 is highest but out of order. Its field contains more GFAP-positive tissue (4.4% of
  voxels vs 2.3–2.8%), which points to field content rather than concentration.
- Projecting the same physical thickness in every stack changed means by ≤ 2%, so unequal
  stack depth is not the explanation.
- The 99.8th percentile is the same for all four (within 1%) and cannot discriminate
  between dilutions (see 4.5).

### 4.2 Signal through the z-stack ("how much darker is the middle?")

- Every stack is **brightest near both tissue surfaces and dimmest at mid-depth**. This fits
  antibody entering a free-floating section from both faces.
- Minimum of the per-slice mean over the middle half of each stack, relative to its
  near-surface peak: **57%, 55%, 63%, 57%** (1:250, 1:500, 1:1000, 1:2000).
- After subtracting background: **43%, 39%, 49%, 41%**.
- Middle-window / first-window ratio of the mean: 0.73, 0.81, 0.87, 0.76.
- **The dip is the same at every dilution.** 1:250, 1:500 and 1:2000 overlap almost exactly
  at mid-depth, so more primary antibody did not improve penetration. The cause of the dip
  (secondary-antibody penetration, optics with a dry objective, tissue structure) remains
  open.

### 4.3 Background fluorescence

| Dilution | Background proxy median | Foreground median | Signal / background |
|---|---|---|---|
| 1:250 | 24 | 542 | 22.6 |
| 1:500 | 21 | 541 | 25.8 |
| 1:1000 | 28 | 558 | 19.9 |
| 1:2000 | 21 | 542 | 25.8 |

- Background is **low (21–28 counts vs ~540 foreground) and does not increase with antibody
  amount**, so there is no sign of concentration-dependent non-specific binding.
- Without no-primary controls, this is a **proxy** (low-intensity tissue away from stained
  structures), not a measurement of autofluorescence.
- It was stable when the threshold was varied 0.5–2× (19–31 counts). The
  signal-to-background ratio, however, depends strongly on the threshold (2.4–2.7-fold
  change), so it is not used to rank dilutions.

### 4.4 Intensity histograms

- Voxel-intensity distributions of all four stacks have **the same shape**: parallel on a
  log scale from 100 to 1,200 counts.
- If more antibody made each structure brighter, the distribution would stretch toward
  higher intensities. It does not.
- 1:1000 is shifted up with the same slope: more GFAP-positive volume, not brighter
  staining.
- No change in the low- vs high-intensity balance across dilutions. Rose's thin/thick-process
  hypothesis is therefore not supported here, although histograms alone cannot measure
  process thickness.

### 4.5 Methodological finding: detector ceiling

- All four stacks show a small pile-up of voxels at 1,300–1,340 counts and an abrupt cut-off
  near 1,420 counts.
- Their maxima agree within 0.8% (1429, 1421, 1419, 1418) although they come from different
  sections.
- This is the signature of **soft saturation (compression) of the photon-counting
  detector** at the current settings. There is no hard clipping: only one voxel sits at each
  maximum.
- **Consequence:** maxima and 99.8th-percentile values measure the detector ceiling, not
  the stain. This includes the p99.8 normalisation and depth-profile method used
  previously, which underestimates the mid-depth dip.
- Comparisons therefore rely on means, medians and the 100–1,200-count range, which about
  99.9% of voxels occupy. The absence of a concentration effect holds there too.

## 5. Overall interpretation and confidence

- **Most likely interpretation:** at 1:2000 this antibody aliquot is already in excess
  (binding plateau) under the current protocol. Diluting from 1:250 to 1:2000 changed none
  of brightness, depth penetration, background or intensity distribution.
- **What can be stated now:** there is no evidence that a primary more concentrated than
  1:2000 helps. 1:2000 gave equivalent results with the least antibody, but this is a
  **tentative observation, not a recommendation**.
- **Confidence: low to moderate.**
  - One field per dilution, from four different sections (field-to-field variation is
    confounded with concentration).
  - Anatomical comparability of the sections is unconfirmed.
  - The brightest signal is compressed by the detector.
  - Background is a proxy only.
  - Depth is nominal (dry objective).
  - The 1:250 stack is about 2 µm shorter.
- The quantitative results are consistent with Rose's visual impression that the four
  dilutions look the same.

## 6. Next steps

1. **Detector linearity test:** image one field at two laser powers or accumulation
   settings. Then choose settings that keep the brightest signal well below about 1,400
   counts.
2. **Replicates:** image several fields per dilution in matched anatomical regions (e.g.
   ipsilateral and contralateral, as Rose planned).
3. **Extend the dilution series** (1:4000, 1:8000) to find where the signal starts to fall.
4. **Include no-primary controls** so background can be measured.
5. **Compare with Emma's acquisition settings**, since her aliquot did show concentration
   dependence.
6. **Analysis additions:**
   - background-subtracted depth ratios as a pipeline output;
   - a morphology check of the deconvolved images against the raw ones;
   - support for multiple fields per concentration.
7. **Compare methods with Rose's earlier Python script** (planned by Rose).
8. **Later project stages (Aims 2–4):** resume napari annotation of Rose's ~50 existing
   tdTomato/GFAP images; segmentation later.

**Questions sent to Rose:**
- Are the four sections from the same region of the cord?
- Is HyD S photon counting expected to be linear at these count levels?
- Did earlier p99.8-based analyses use similar acquisition settings?

## 7. Figures available (PNG, in the repository under `reports/figures/`)

1. **Max projections of the four dilutions on one shared display scale.** They look alike.
   1:1000 is densely stained across the whole field; the others are sparser in their lower
   half.
2. **Per-slice mean vs depth.** All four curves peak near 3–4 µm, dip at 11–13 µm and rise
   again toward the second surface. 1:250, 1:500 and 1:2000 nearly overlap; 1:1000 is higher
   throughout.
3. **Per-slice 99.8th percentile vs depth.** Rose's previous method. Surface peaks are
   flattened at about 1,300 counts by the detector ceiling.
4. **Threshold-sensitivity diagnostic.** Background is stable across thresholds;
   signal-to-background is not.
5. **Mask QA.** Foreground (shared threshold) and background-proxy regions overlaid on the
   max projections.
6. **Voxel-intensity histograms.** Identical bins, log scale. Parallel curves, with a
   pile-up and cut-off near 1,400 counts in all four.
7. **Survival curves** (fraction of voxels at or above each intensity). 1:250, 1:500 and
   1:2000 nearly coincide; 1:1000 is shifted up; all converge at the ceiling.

## 8. Deliverables produced

- `gfap-conc` analysis pipeline: Python, 90 tests, documented, in a private GitHub
  repository.
- English results report for Rose (Markdown and .docx, 9 pages, 7 figures):
  `reports/GFAP_concentration_curve_report_2026-11-23.*`.
- Full output snapshot of the confirmed run (all tables, figures, QC warnings, provenance
  record): `reports/run_outputs/run_20261122-160023/`.
- Methods, validation status, session log and an AI-handoff document in `docs/`.
