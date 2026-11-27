# AI handoff: GFAP concentration-curve analysis

**Purpose of this file:** give a new AI assistant everything needed to continue this
project without the original chat. Read it fully before acting. Last updated 2026-11-23.

## 1. TL;DR

- **Project.** Quantify a Mouse anti-GFAP primary-antibody concentration curve (1:250,
  1:500, 1:1000, 1:2000) from Leica confocal z-stacks. This is a PHAR 599 student project
  (McGill, Ribeiro-da-Silva lab). Student: **Sewon** (the user). Mentor: **Rose Rodrigues**.
  PI: Prof. Alfredo Ribeiro-da-Silva.
- **Question.** Rose saw identical brightness at all four dilutions (unlike Emma's aliquot)
  and asked for four analyses: (1) average bit value of max projections, with a 99.8th-
  percentile "normalization max"; (2) intensity change through the stack; (3) background
  fluorescence; (4) intensity histograms.
- **Status.** A Python pipeline (`gfap-conc`) is built, tested (92 tests) and has been run on
  the real file. A first report exists. **Result:** no concentration dependence in any
  measure. There is also soft detector saturation near 1,400 photon counts.
- **Next step.** Waiting for Rose's feedback. Likely next work: background-subtracted
  depth ratios, a linearity analysis for re-imaging, `_Lng` morphology QA, and support for
  several fields per concentration (section 9).

## 2. Read in this order

1. This file.
2. `docs/context/03_Coding_Task_and_Analysis_Spec.md`: the original spec. Rose's emails are
   quoted verbatim in all three context files.
3. `reports/GFAP_concentration_curve_report_2026-11-23.md`: the results as sent to Rose,
   with figures in `reports/figures/`. **The complete raw output of the confirmed run**
   (all CSVs, 12 figures, `summary.md`, `qc_warnings.txt`, `run_manifest.json`,
   `config_used.yaml`) is in `reports/run_outputs/run_20261122-160023/`. Use it to check or
   extend results instead of asking for re-runs. The files were written on Windows (CRLF;
   `config_used.yaml` has a UTF-8 BOM).
4. `docs/methods.md`: exact definitions of every number.
5. `docs/VALIDATION_STATUS.md`: what is tested, and what was learned from the real data.
6. `docs/SESSION_LOG.md`: the history, every problem fixed, and an authorship note.
7. `README.md` and `docs/RUNNING_ON_WINDOWS.md`: how to run.

## 3. Rules that must not be broken

- **Never modify the raw `.lif`/`.lifext`.** Open them read-only. Every run checks the SHA-256.
- **Quantify raw images only.** `_Lng` images (deconvolved, rescaled to the full 16-bit
  range) are blocked from all intensity analyses by `require_raw()`. Rose: they are not for
  answering "which concentration works best".
- **Do not rescale before comparing.** p99.8 is a display ceiling and a separate column.
  Per-image normalised values are secondary only.
- **Use one shared absolute foreground threshold** for all concentrations, never
  per-image thresholds.
- **Background is a proxy unless controls exist.** Never call it autofluorescence.
- **Review is mandatory** before comparing concentrations. `gfap-conc run` refuses unless
  `review:` in the config is confirmed with the hash of the current mask and z-range
  settings.
- **Statistics are descriptive.** There is one field per concentration, so no p-values;
  pixels and slices are not replicates. "No defensible optimum yet" is an acceptable
  conclusion.
- **Keep synthetic and real outputs apart:** `results/synthetic/` vs `results/real/`, and
  the `data_origin` column.
- **Do not invent results.** If data are not visible to you, ask the user to run a command
  and paste the output.
- **Authorship.** The code was written by an AI (Claude) under Sewon's direction. Rose
  originally wanted Sewon's own implementation to compare with her script, so this must be
  disclosed. It already is, in the report and in `SESSION_LOG.md`. Keep it that way.

## 4. Environment and how to work with the user

- **Code lives in the GitHub repo** `sjung-research/phar599` (private), branch
  `claude/blissful-cori-wjn019`.
- **The user runs everything on a Windows 11 laptop in PowerShell:**
  - repo at `C:\Users\sjung3\Documents\phar599`;
  - data at `E:\2026NovemberGFAPlif\`;
  - Korean locale (cp949).
- **A second Windows computer (the lab desktop) also works.** There the data are at
  `C:\Users\sjung3\Desktop\OneDrive_1_11-23-2026\` (same SHA-256). A run there on
  2026-11-23 reproduced all nine CSV tables byte-for-byte, as checked by the student.
- **The data never leave the laptop.** A cloud AI only sees what the user pastes: console
  output, CSV excerpts, `summary.md`, PNG figures.
- **Workflow:** the AI edits code → pushes → the user runs `git pull`, then
  `uv run pytest`, then the command.
- **User preferences:**
  - explanations in **Korean**;
  - program output, reports and docs in **English**;
  - give **one command per code block**: the user once pasted two commands on one line,
    and once pasted YAML into PowerShell;
  - prefer commands that write files (PowerShell here-strings, `-replace`) over asking the
    user to edit YAML by hand.
- **Windows gotchas already handled:**
  - all file I/O is explicit UTF-8;
  - config loading tolerates a BOM;
  - console output never crashes on non-ASCII characters;
  - Windows paths in YAML go in **single quotes**.
- **Tools:** uv; Python 3.13 (`.python-version`); dependencies pinned in `uv.lock`.

## 5. Data facts (confirmed from metadata)

| Item | Value |
|---|---|
| File | `E:\2026NovemberGFAPlif\Sewon_GFAP_curve.lif`, 2,403,467,585 bytes, SHA-256 `da134597fa839257395a3f5595f4c0a634910837337f518725bc278fc27041be` |
| Companion | `Sewon_GFAP_curve.lifext`, 384,575,069 bytes, SHA-256 `534d60073c64663a78426e7efa98944217e695dc56be60deb896d9732db5e1a3`; contributes no images through the `.lif`; not inspected separately yet |
| Series | `1-250`, `1-500`, `1-1000`, `1-2000` (raw, ZCYX) and `1-250_Lng` … `1-2000_Lng` (CZYX, max 65535). **The suffix is `_Lng`, not Rose's `_lng`.** |
| Channels | 0 = DAPI (405 nm, 2%), **1 = ALEXA 647 (GFAP)** (638 nm, 2.5%) |
| Geometry | 2200 × 1500 px at 0.13216 µm; z-step 1.0007 µm; 21 / 24 / 23 / 23 slices |
| Acquisition (identical for all four) | HC PL APO CS2 20x/0.75 **dry**, zoom 2; HyD S (SiPM) **PhotonCounting**, gain 2.5, offset 0; scan speed 400; line accumulation 4; frame average 2; pinhole 1.0 AU; `Resolution` 16 bit |
| Differences | Stage X 66.6 / 62.2 / 57.8 / 53.2 mm (four different sections); z start/end; slice count |
| Value range | Raw AF647 max 1429 / 1421 / 1419 / 1418; one voxel at each max; background median about 20 counts |

## 6. The exact config used for the confirmed run

`config/analysis.yaml` is git-ignored. This is its content; an exact copy is
`reports/run_outputs/run_20261122-160023/config_used.yaml`. The `analysis:` block was
omitted, so all defaults apply. Defaults at the time of the run were as in
`config/analysis.example.yaml`, except two later changes that do not affect the review hash:
`histogram.range_max: observed` and `qc_max_fraction_of_tissue: 0.97`.

```yaml
input:
  lif_path: 'E:\2026NovemberGFAPlif\Sewon_GFAP_curve.lif'
  lif_sha256: da134597fa839257395a3f5595f4c0a634910837337f518725bc278fc27041be
  gfap_channel: 1
  reference_channel: 0
  bit_depth_override: null
  intensity_unit: 'photon counts (HyD S PhotonCounting; line acc. 4, frame avg. 2)'

series:
  - {name: '1-250',      concentration: '1:250',  kind: raw}
  - {name: '1-500',      concentration: '1:500',  kind: raw}
  - {name: '1-1000',     concentration: '1:1000', kind: raw}
  - {name: '1-2000',     concentration: '1:2000', kind: raw}
  - {name: '1-250_Lng',  concentration: '1:250',  kind: deconvolved, pair_of: '1-250'}
  - {name: '1-500_Lng',  concentration: '1:500',  kind: deconvolved, pair_of: '1-500'}
  - {name: '1-1000_Lng', concentration: '1:1000', kind: deconvolved, pair_of: '1-1000'}
  - {name: '1-2000_Lng', concentration: '1:2000', kind: deconvolved, pair_of: '1-2000'}

review:
  tissue_mask_confirmed: true
  z_range_confirmed: true
  reviewed_settings_sha256: 6fd16a9c650c08a97525ff47892088c32e02a03f9756c235af9e2fdc7bd9b3cb
  reviewed_by: Sewon
  notes: 'Mask none: all four fields fully inside tissue. z-range full: no empty slices; darker end slices (1:500 z22-23, 1:1000 z22) are tissue surfaces. 1:250 stack spans 20.0 um vs 22.0-22.3 um and likely stops short of the lower surface, so middle_over_last is not strictly comparable; use middle_over_first. Four different sections (stage X 66.6/62.2/57.8/53.2 mm); region comparability to be confirmed with Rose.'

output_dir: ../results/real
```

The review decisions were tissue mask `none` (every field lies fully inside tissue) and
z-range `full` (no empty slices). They were made on `review_tissue_mask.png` and
`review_z_range.png`. Changing any mask or z-range setting invalidates the hash: re-run
`gfap-conc review` and re-confirm.

## 7. Pipeline in one screen

Workflow commands (PowerShell):

```text
uv run gfap-conc inspect "<lif>"                  # inventory, metadata, checksum, config suggestion (read-only)
uv run gfap-conc review --config config\analysis.yaml   # mask/z-range QA figures + checklist; prints settings hash
uv run gfap-conc run --config config\analysis.yaml      # the four analyses; refuses without confirmed review
uv run gfap-conc demo                             # synthetic end-to-end run (watermarked)
```

**Run outputs** go to `results/real/run_<timestamp>/`:
- `max_projection_summary.csv`, `depth_profile.csv`, `middle_vs_outer_summary.csv`;
- `background_summary.csv` (with `bg_qc_flags`), `threshold_sensitivity.csv`;
- `intensity_histogram.csv`, `intensity_survival.csv`, `histogram_quantiles.csv`;
- `qc_stack_summary.csv`, `qc_warnings.txt`;
- `figures/*.png`, `summary.md`, `run_manifest.json`, `config_used.yaml`.

**Code map** (`src/gfap_conc/`):
- `model.py`: `GFAPStack`, `require_raw`.
- `stats.py`.
- `projection.py`: analysis 1.
- `depth.py`: analysis 2 and the auto z-range rule.
- `masks.py`.
- `background.py`: analysis 3, QC flags, sensitivity.
- `histograms.py`: analysis 4.
- `pipeline.py`: orchestration and QC warnings.
- `io_lif.py`, `loading.py`: `.lif` access.
- `review.py`.
- `plots.py`.
- `report.py`.
- `config.py`: YAML schema; unknown keys are rejected.
- `cli.py`.

Tests are in `tests/`. `tests/lif_writer.py` writes minimal LIF files.

**Report:** `reports/GFAP_concentration_curve_report_2026-11-23.md` is the source. Build the
docx with:

```text
node reports/build_report_docx.js reports/<name>.md reports/<name>.docx
```

There is also an online copy (Claude Docs): <https://claude.ai/artifact/EMni7yBiV8BouX5G2bJP55>.

## 8. Results (confirmed run `run_20261122-160023`, code `e5552bf`, photon counts)

| Dilution | MIP mean | MIP p99.8 | Stack median / q0.9 / q0.99 | Mid/first window (mean) | Background proxy | SBR @ T=341 | GFAP voxels ≥ T |
|---|---|---|---|---|---|---|---|
| 1:250 | 219.5 | 1357 | 27 / 126 / 662 | 0.73 | 24 | 22.6 | 2.8% |
| 1:500 | 204.6 | 1352 | 22 / 106 / 590 | 0.81 | 21 | 25.8 | 2.3% |
| 1:1000 | 288.3 | 1364 | 33 / 173 / 849 | 0.87 | 28 | 19.9 | 4.4% |
| 1:2000 | 199.4 | 1354 | 23 / 112 / 609 | 0.76 | 21 | 25.8 | 2.4% |

**Interpretation, as in the report:**
- There is no dose dependence. The antibody is likely in excess (plateau) at ≥ 1:2000.
- 1:1000 differs in GFAP-positive *volume* (field content), not in brightness: its
  histogram has the same slope, shifted up.
- The mid-depth dip is the same at every dilution. The minimum of the per-slice mean over
  the middle half of each stack is 57 / 55 / 63 / 57% of the near-surface peak. After
  subtracting the per-slice `bg_proxy_median`, it is 43 / 39 / 49 / 41%. These were computed
  ad hoc from `depth_profile.csv` and are not yet a pipeline output.
- **Soft photon-counting saturation:** there is a bump at 1,300–1,340 counts and a cliff at
  about 1,420 in every stack. Maxima and p99.8 measure the ceiling, so comparisons use
  means, medians and the 100–1,200-count range.
- Confidence is low to moderate: n = 1 field per dilution, four different sections, and
  region comparability is unconfirmed.

## 9. Open tasks (suggested order)

1. **Wait for Rose's answers.** Are the sections the same region? Is HyD S photon counting
   linear at these counts? Were her earlier p99.8 analyses acquired the same way?
2. **Background-subtracted depth retention**, e.g. `(metric − bg_proxy_median)` ratios
   using the per-slice `bg_proxy_median` already in `depth_profile.csv`. Implement in
   `depth.py` and `pipeline.py` as additional columns, keeping the raw ones. The report's
   43 / 39 / 49 / 41% values were computed ad hoc (peak = max of the first half,
   trough = min over the middle 50% of slices). Make this a tested pipeline output.
3. **Linearity analysis for re-imaging:** if Rose images one field at two laser powers or
   accumulations, compare them pixel by pixel (ratio vs intensity) to locate the onset of
   compression. Needs a new small module and a CLI command.
4. **Several fields per concentration** (ipsi/contra, replicates). The pipeline currently
   refuses duplicate concentrations on purpose. It needs an explicit aggregation design:
   per-field results, then a summary across fields.
5. **Done 2026-11-27:** figures in the layout of Rose's first-staining PDFs
   (`mip_with_histograms_rose_style.png`, `depth_slice_gallery_{per_image,shared}.png`;
   `plots.py`, `depth.gallery_slice_indices`). They are not in the committed run snapshot
   until the analysis is re-run. **Finding from comparing with Rose's PDFs:** her Mouse and
   Chicken MIP histograms show the same pile-up and cliff (about 1,000 and 1,050 counts, at
   63x), with p99.8 1005–1018 for all six images. This independently supports the
   detector-ceiling interpretation. Rabbit, which is dimmer, shows no pile-up.
6. **`_Lng` morphology QA:** paired raw vs `_Lng` crops, each on its own display range. The
   config already has `pair_of`.
7. **Inspect the `.lifext` separately**
   (`uv run gfap-conc inspect "E:\...\Sewon_GFAP_curve.lifext"`) to know what it holds.
8. **Later, out of scope now:** napari annotation of Rose's ~50 earlier astrocyte images, and
   segmentation.

## 10. Pitfalls already met (do not repeat)

- **Histogram bins over the nominal 16-bit range** are 256 counts wide. Default is now the
  observed range.
- **The original auto z-range rule** (half-max of the profile) cut the real stacks at their
  mid-depth dip. It now uses the background floor and only trims the ends. Still prefer
  `manual` ranges chosen in review.
- **The `LARGE_REGION` flag at 85%** fired for all stacks: GFAP is sparse at 20×. It is now
  97%.
- **Per-slice p99.8 depth curves** (Rose's earlier method) are flattened by the detector
  ceiling at the surfaces.
- **Do not assume Rose's naming:** `_Lng` vs `_lng`; dye names are in the metadata.
