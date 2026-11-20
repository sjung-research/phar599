# 03 — Claude 코딩 작업 명세: Mouse-GFAP concentration-curve `.lif` 정량 분석

> **우선순위 가장 높은 문서.** 이 문서를 Claude에게 첨부하고 `.lif` 실제 원본을 전달해 Python code, reproducible report, plots를 요청한다. 역사적 배경은 `01_Background_and_Methods.md`, 큰 연구 목적은 `02_Current_Research_Project.md` 참조.
>
> **2026-11-19 Rose의 요청에서 확실한 것:** (1) max-projected image average bit value + 99.8th-percentile normalization maximum, (2) z-stack 내 intensity 변화, (3) background fluorescence 차이, (4) intensity histogram 분포, (5) `_lng` images는 deconvolved + full 16-bit range normalized → 농도별 raw intensity 비교에 쓰지 말 것, (6) 학생의 Python code와 Rose의 이전 Python code를 나중에 비교한다.
>
> **현 상태:** Rose가 SharePoint 접근 이메일을 보내겠다고 했다. **실제 `.lif` 파일 내용/이름/크기/메타데이터는 아직 이 문서 작성 시점에 확인되지 않았다.** 아래에는 구현 설계 및 QC 조건만 있다. 불명확한 것을 확정된 사실처럼 hardcode하지 말 것.

## 0. Claude에게 줄 직접 명령 (task statement)

**목표:** 원본 Leica `.lif` fluorescence z-stack에서 Mouse anti-GFAP primary concentrations `1:250`, `1:500`, `1:1000`, `1:2000`를 비교하여 **어떤 농도가 구조, 신호 강도, 깊이 균일성, background를 균형 있게 만족하는지** 정량 평가하는 Python pipeline을 설계·구현한다.

**필수 산출물:** 재실행 가능한 Python source code, 환경/사용법 문서, machine-readable numeric CSV, 농도별 비교 figure들, 결과 해석과 한계가 명시된 brief report. **원본 `.lif`는 수정하거나 덮어쓰지 말 것.** 반드시 `raw` 원본 image series를 우선 분석한다. `_lng` series는 morphology cross-check용 보조 패널로만 사용한다. Rose의 과거 code를 보지 않은 상태에서 먼저 독립적인 분석안을 만들고, 나중에 code comparison이 가능하게 한다.

**진행 원칙:** 자료를 못 읽으면 가짜 값을 출력하지 말고 import 단계/미확인 metadata를 보고한다. 원본이 도착하면 먼저 series, 채널, z축, dtype, acquisition metadata를 열람한 뒤 실제 파일 구조에 맞게 코드를 작성한다. 데이터 구조와 계산 정의를 보고서에 명시한다.

---

## 1. 실험 design과 핵심 비교

### Experimental conditions

| Field | Value / notes |
|---|---|
| Primary | Mouse anti-GFAP (newly acquired antibody corresponding to Rose's previous successful product) |
| Independent variable | Primary dilution: `1:250`, `1:500`, `1:1000`, `1:2000` (2-fold serial dilution) |
| Secondary | anti-Mouse AF647; lab communication confirms AF647 + DAPI for Day 2 |
| Sample | mouse spinal cord fixed tissue, typically ~50 µm physical sections; recorded z span may differ |
| Source data | Leica `.lif` image container from SharePoint |
| Images in `.lif` | Rose: one image per concentration, plus another image per concentration with suffix `_lng` |
| `_lng` | deconvolved + normalized to entire 16-bit range; **DO NOT use for identifying best concentration by absolute intensity** |
| Raw replicates | Rose only explicitly promises **one raw image per concentration** in this `.lif`; do not infer n=3 or biological replicates |
| ROI / side | Rose previously contemplated ipsilateral and contralateral images; check whether delivered `.lif` actually contains those |
| Imaging parameters | Need to inspect `.lif` metadata; equal detector/laser/scale settings must be assessed, not assumed |

### Real scientific question

The staining worked, but Rose saw approximately **the same apparent brightness at all four concentrations** for this aliquot; Emma's sample varied with concentration. We need to test quantitatively whether (a) there really are minimal differences, (b) brightness might be similar but depth/background/histogram differs, or (c) acquisition/display/ROI characteristics obscure a concentration trend. A rigorous outcome may be **“no defensible optimal concentration from currently available images.”**

### Distinguish research questions

- **NOW:** Quantitative intensity/penetration/background/histogram analysis of raw concentration-curve images.
- **LATER:** Annotate pre-existing ~50 astrocyte images in napari; eventual segmentation model.
- **NOT requested now:** Train a GFAP segmentation model or compare biological treatment groups. Do not allow the project background to expand the scope unnecessarily.

---

## 2. Workflow stage A — inspect `.lif` and preserve provenance

Before scientific calculations:

1. **Inventory** the `.lif` image series: exact series names, `raw` or `_lng`, candidate concentration, number of channels, number of z-planes, dimensions, dtype, valid range, physical z-step (`µm`), XY pixel sizes, microscope/scanner metadata, acquisition timestamps if available.
2. **Choose a capable `.lif` reader:** e.g. Bio-Formats-compatible reader or appropriate Python LIF reader. Prefer one that exposes OME-compatible axis metadata. Check software compatibility rather than assuming arbitrary `.lif` files have the same layout.
3. Verify channel identities: GFAP **AF647** channel (and possibly DAPI channel). **Do not quantify GFAP from DAPI** or from a merged RGB preview.
4. Convert axis order to explicit labelled `(Z, Y, X)` for the selected GFAP channel; assert dimensions before statistics. `T` and `C` should be handled explicitly if present.
5. Pair original with `_lng` without altering either, retaining the original names and metadata.
6. Inspect raw intensity limits and saturation fraction (e.g. how often exactly full-scale maximum appears). Record what constitutes true bit depth versus displayed `vmin/vmax`.
7. Preserve **immutable raw data**; results go to a separate output directory.

**Proposed initial output:** `series_inventory.csv` with at least `series_name, concentration, image_kind, channel, shape, Z, Y, X, dtype, z_step_um, z_extent_um, XY_pixel_size, raw_min, raw_max, saturation_fraction, acquisition_settings, provenance_notes`.

**Important:** If the 4 images have different z-depth, XY field, image size, detector gains, laser power, resolution, or section/animal region, report it visibly. Some numeric comparisons may then need restriction to matched ROIs or acquisition-matched subsets.

---

## 3. Workflow stage B — four required analyses from Rose

### Analysis 1 — max projection: **average bit value** + 99.8th-percentile display reference

For each raw GFAP stack `(Z,Y,X)`:

- Compute **MIP**: `MIP[y,x] = max_z I[z,y,x]`.
- Compute **RAW average intensity**: arithmetic mean of the raw MIP pixels, **reporting numeric units from source** (digital counts/ADU unless calibrated otherwise).
- Compute `median`, `std`, `p95`, **`p99.8`**, `min/max`, fraction 0-valued and fraction saturated. Optional: masked/ROI-restricted versions.
- Make a **magma-like** grayscale/color figure, with a labeled display range; choose Rose-like `vmax = p99.8(raw MIP)` for a **within-image display view**, and add a second version with a **shared display scale across concentrations** for honest visual comparison.
- Keep `raw_mean` unscaled, separately report `p99.8` as **display normalization max**; optionally provide normalized mean as a *secondary* relative shape/contrast diagnostic, not the only cross-concentration outcome.
- **Warning:** per-image normalization can make a weak image look as bright as a strong one. Never use individually rescaled image values alone to decide which dilution is brighter.
- **Warning:** MIP mean can depend on `Z`/depth coverage; include a per-plane or matched-depth analysis alongside it.

**Proposed outputs:** `max_projection_summary.csv`; `figures/raw_mip_grid_shared_scale.png`, `figures/raw_mip_grid_p998_display.png`, optional TIFF outputs preserving source intensities, figure captions documenting display transforms.

### Analysis 2 — **intensity changes across the z-stack** (is the middle darker?)

Per raw GFAP z-slice, calculate and tabulate:

- depth index `z`, physical depth `µm` from metadata, normalized depth fraction `z/(Z-1)`;
- `mean`, `median`, `p95`, `p99.8`, `std`, **spatial pixel variance within the slice** (if helpful), background-adjusted stats **only if a defensible background estimate is available**;
- optional foreground fraction above a consistent, pre-defined signal threshold. Do not make each concentration have its own adaptive foreground threshold without explaining the consequence.

Compare:

- plot `p99.8` versus physical `z` and/or normalized depth (Rose's previous GFAP depth-profile figure style);
- plot **raw mean intensity per slice** versus z;
- report center-versus-edge (or middle-versus-outer) signal ratio, e.g. middle 20% z slices versus first/last 20% under a predeclared definition; choose thickness-aware windows rather than assuming all stacks are identical;
- optional depth CV and range, and top-to-center decline. Note that one edge can be brighter because of orientation/geometry; do not assume symmetric profiles.

**Critical distinction:** Variation of a depth profile is **not the same** as within-slice pixel variance. Rose previously said *“variance over stack”* in a broad sense, and on Nov 19 specified depth intensity changes; one can report both clearly defined quantities, never equate them.

**Proposed outputs:** `depth_profile.csv` (`concentration,z,z_um,mean,median,p99_8,variance,...`), `depth_profile_raw_p998.png`, `depth_profile_raw_mean.png`, `middle_vs_surface_summary.csv`.

### Analysis 3 — background fluorescence differences

Rose says bad perfusion can cause tissue autofluorescence; **the actual cause is unverified**, and not all low-intensity pixels are background.

- Ask/inspect if there are **no-primary controls** or annotated negative ROIs in the delivered data.
- If good background regions can be identified (same tissue/section orientation & tissue ROI strategy), calculate **background median, mean, robust dispersion**, and optional signal-to-background or contrast metrics.
- If no controls or convincing ROIs exist: provide exploratory **low-percentile/low-signal-region proxies** clearly marked `background proxy`, not a definitive autofluorescence measurement.
- Compare background in standardized units; do **not** manipulate individual series to force equal backgrounds.
- Include ROI/mask visual QA images; avoid treating all black pixels *outside* the specimen as biological tissue background.

**Proposed outputs:** `background_summary.csv`, `figures/background_ROI_overlays.png`, `figures/background_comparison.png`, optional background-corrected summaries if scientifically justified.

### Analysis 4 — pixel-intensity histogram distribution

- Build raw grayscale GFAP intensity histograms for each concentration, preferably with **identical bin edges and range** for direct comparison.
- For 16-bit source, use configurable bin counts (e.g. 256/512/1024 depending on distribution); perhaps zoom into non-saturated lower-intensity region and compare with full-range panel.
- Show normalized frequency (or density) and optionally log-count y-axis; avoid hiding near-zero counts via display choices.
- Provide MIP histograms and, if practical, 3D-stack aggregate pixel histograms **as separate outputs** (they answer different questions).
- If reasonable tissue/foreground masks are available, separate all-pixel histogram from tissue ROI histogram.
- Rose hypothesized: many **thin processes** may correspond to more low-intensity counts, thicker/brighter processes to higher-intensity distribution. **Histogram shape alone cannot establish process thickness**; discuss this as a morphological hypothesis rather than proven cell-process classification.

**Proposed outputs:** `intensity_histogram.csv`, `figures/histograms_raw_shared_bins.png`, `figures/histograms_zoomed.png`.

---

## 4. `_lng` deconvolved images — separate, morphology-only comparison

The suffix is `_lng` **exactly as Rose typed**; do not assume it stands for a specific proprietary Leica algorithm or normalize it to another suffix.

1. Match each `_lng` image with original by name and plausible metadata/field-of-view, record mapping explicitly.
2. Show matched visual crops / process visibility side-by-side, using clearly noted **separate display ranges** if needed.
3. Optionally compare structural similarity after explicitly stated registration/normalization, but do not treat pixel equality as expected.
4. **Do not compute “best concentration” using normalized `_lng` raw histogram mean, p99.8, background, or absolute brightness.** Their intensities have been rescaled to occupy the full 16-bit range.
5. Report any discrepancies (e.g. new structures, oversharpening, ringing, apparent changes) cautiously without implying an independent new raw photon measurement.
6. Rose said later the lab may work **exclusively with deconvolved versions** once the question is different (segmentation morphology). That future plan does **not** override the current concentration-curve restriction.

---

## 5. Recommended Python project structure and deliverables

```text
GFAP_concentration_analysis/
├── README.md                      # setup, data location, run commands, interpretation
├── requirements.txt / environment.yml
├── config.yaml                    # input paths, GFAP channel, series mappings, options
├── src/
│   ├── inspect_lif.py             # list series, axes, metadata, source-image QC
│   ├── image_io.py                # raw/deconvolved loading and shape validation
│   ├── quantification.py          # MIP, depth profile, background, histograms
│   ├── visualize.py               # consistent figures, legends and color scales
│   └── analysis.py                # reproducible run entrypoint
├── tests/                          # small synthetic stacks, no lab data required
│   └── test_quantification.py
└── results/
    ├── series_inventory.csv
    ├── max_projection_summary.csv
    ├── depth_profile.csv
    ├── middle_vs_surface_summary.csv
    ├── intensity_histogram.csv
    ├── background_summary.csv      # only valid/labelled background estimates
    ├── figures/
    └── summary.md                   # findings, uncertainty, QC notes
```

This tree is a **suggested implementation**, not Rose's required file naming. Claude may change library/module structure if it improves maintainability.

### Minimal tests before analyzing real data

- A synthetic 3D stack with known slice means, `p99.8`, variance produces expected values.
- MIP equals elementwise max over the **correct z axis**.
- Series containing `T/C/Z` are handled explicitly; no unnoticed channel mixing.
- A bright single outlier changes `max` far more than `p99.8`, as expected.
- Identical distributions at different bit-depth scales are **not mistakenly considered equal absolute intensity** after per-series normalization.
- `_lng` names are excluded from **primary concentration comparison**; included only in morphology QA.
- Histogram bins are shared across conditions; CSV and plotted totals agree within numerical tolerance.
- Background ROI coordinates align with original XY and are not accidentally applied after resizing without mapping.
- Output retains association between concentration, actual series name, GFAP channel, and every statistic.

### Output interpretation template

For each concentration, explain:

- `Raw MIP mean` and `p99.8` (units, display assumptions)
- `Center vs outer-depth signal retention`
- `Background / signal-to-background` (if reliably measurable)
- `Histogram shape and plausible morphological interpretation`
- `Acquisition/dtype/ROI comparability`, sample count
- `Conclusion` and confidence: strongest evidence, trade-offs, unknowns

A concentration with the largest mean intensity is **not automatically** the best: high background, poor depth retention, unresolved/overlapping morphology, saturation could outweigh brightness. With one image per condition, **avoid unqualified p-values, pseudo-replication, or statements that the concentration causes specific effects**. Technical image pixels/z-slices are not independent biological replicates.

---

## 6. Ambiguities Claude must not silently resolve

| Unknown | Required handling |
|---|---|
| Actual content of `.lif` from SharePoint | Inspect and inventory, report if missing |
| Exact raw/deconvolved names and concentration mapping | Use metadata / names; seek confirmation when ambiguous |
| GFAP channel, axis order, `z` step, bit depth | Read metadata and verify by inspection |
| Absolute intensity comparable across acquisitions? | Compare gain/laser/detector/ROI settings; qualify if not |
| Rose's exact normalization formula | Her `normalization max` likely p99.8 used for display; keep raw means separate. Ask if necessary |
| What “variance over stack” should denote | At least explicit depth variation; optionally separate per-slice pixel variance and variance of depth-series summary |
| Background region/control | Label actual vs proxy. Do not call a low percentile autofluorescence without evidence |
| How many true biological replicates? | She states one raw series per concentration in this `.lif`; no inferential statistics unless more exist |
| Ipsi/contra additional imaging | Only analyze separately if source files actually contain identifiable side metadata |
| Right toolchain / `.lif` reader package versions | Confirm compatibility with actual file; do not invent a successfully executed run |

---

## 7. Reference figures and why they matter

The user may separately attach:

- `GFAP_maxproj_all_magma.pdf` — a prior **Mouse/Rabbit/Chicken** comparison with representative max projections and histogram plots, three triplicates per antibody. **Not** current concentration curve.
- `GFAP_depth_slices_magma.pdf` — prior stack depth profiles using **per-slice 99.8th percentile** and representative z-slice gallery. **Not** current concentration curve.

Use them as a **figure-output style guide** and conceptual baseline. Do not copy their measured numeric results into the new Mouse-only concentration analysis, and do not assume Rose's Python source is bundled with the PDF.

---

## 8. Rose 이메일 원문 A — independent Python assignment (2026-11-18, screenshot)

```text
I'm losing the fight with conda, and I need to finish packing my house for my move Saturday so I have to tap out for this week. In the interim since I think Friday is no-go, how about I give you access to the .lif file for the images that I took of the concentration curve, and you try creating/running the analysis of the bit values and variance over stack yourself?

Are you okay with playing around with that this week while I'm moving houses and getting set up? We can let you try it on your own first, then we can compare python files of how I did it last time vs what you make on your own!

Let me know what you think!
```

## 9. Rose 이메일 원문 B — EXACT scientific analysis directions (2026-11-19, screenshot)

> **가장 중요한 원문.** 아래 영어는 화면에서 읽히는 본문을 가능한 한 그대로 전사했다. 이메일 화면을 읽어 옮긴 것이므로 극소수의 공백/인용부호 glyph 차이는 있을 수 있지만, 요청 항목·표현·주의사항은 원문대로 보존했다.

```text
So, from the images, there’s a few things you can ‘analyze’:

- Average bit value on max-projected images (and I add a ’normalization max’ to 99.8th percentile, because one single crazy high pixel can skew data)
- Intensity changes over the course of the stack (how much darker is it in the middle?)
- Differences in background fluorescence (bad perfusion can cause auto-fluorescence of tissue)
- Histogram distribution of image intensity (if one image has a lot more ’thin’ processes, you’ll have a higher peak at lower intensity fluorescence, thick processes skew towards the top range of the histogram)

The .lif file has one image of each concentration, then another version of the image with _lng on the name - these are deconvolved versions of the images, which have higher resolution but they are normalized to the entire 16-bit range, so those are not what we want to use to answer the question of ‘which antibody concentration works best’. I use them to discern how close they look to the original images (down the line we’ll work exclusively with the deconvolved versions).

I’ll send over an email from Sharepoint with access to the .lif file!
```

## 10. Rose 이메일 원문 C — why the concentration question remains open (2026-11-13)

```text
Hi Sewon,

I took a look at your slides, and somehow, your aliquot of mouse anti GFAP looks exactly the same in all concentrations (unlike Emma’s, which changes intensity based on concentration). I want to go back and take more images from both ipsi and contra sides, to get actual quantification. So I don’t have a “final concentration” for you yet, though the stain absolutely worked.

BUT! We’re just about ready to start staining new samples. I just need to 1) do tamoxifen injections and 2) collect tissue. This will take a couple of weeks minimum, probably more because I’m currently in the middle of packing up to move apartments and poor Tesla still has more upcoming vet appointments.

I think in the meantime, since we’re keeping the same primary antibody as the original experiment, we can actually get into labelling the data from my original 50 images!!!!!!!!!!

On my calendar, I put in next Tuesday to (attempt to) set up an exact replica of my Napari environment on the lab computer, so that I can teach you how to open the image files and use the hotkeys I created in order to fix the labels on the images. If I’m successful with that Tuesday, I could do Wednesday or Friday to show you how it all works. How does that sound?

Sorry for the long email!
```

---

## 11. Final deliverable checklist for Claude

- [ ] 실제 `.lif` 파일을 읽고 series/raw/deconvolved 목록을 **검증**했다.
- [ ] raw Mouse-GFAP concentration series 각각을 올바른 GFAP 채널로 추출했다.
- [ ] `_lng` image를 raw-intensity comparison에서 배제했다.
- [ ] raw MIP mean과 p99.8 display max를 **분리**하여 계산했다.
- [ ] depth profile과 middle-depth attenuation을 비교했다.
- [ ] background (검증된 ROI 또는 명시적 proxy) 비교를 수행했다.
- [ ] 동일 bin을 사용하는 intensity histogram을 만들었다.
- [ ] 모든 figure에 raw/deconv, scale, intensity units, depth units, concentration을 표시했다.
- [ ] 원본을 덮어쓰지 않았고 CSV/figures/summary와 재실행 코드를 제공했다.
- [ ] 분석상 모호함·sample limitations·potential confounds를 기록했다.
- [ ] 결론을 무리하게 특정 농도로 만들지 않고 Rose와 검토할 다음 질문을 정리했다.
