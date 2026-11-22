# GFAP concentration-curve analysis (`gfap-conc`)

Quantitative analysis of a **Mouse anti-GFAP primary-antibody concentration curve**
(1:250, 1:500, 1:1000, 1:2000; Donkey anti-Mouse AF647 secondary) imaged as Leica
`.lif` confocal z-stacks. PHAR 599, Ribeiro-da-Silva lab (McGill), with Rose Rodrigues.

The scientific question is why staining looked equally bright at all four dilutions,
and which dilution best balances signal, depth penetration, background and
morphological detail. The pipeline implements Rose's four requested analyses on the
**raw** (non-deconvolved) images:

| # | Analysis | Main outputs |
|---|---|---|
| 1 | Max-intensity projection: raw mean ("average bit value") + 99.8th-percentile display reference | `max_projection_summary.csv`, `figures/raw_mip_grid_*.png` |
| 2 | Intensity through the z-stack ("how much darker is the middle?") | `depth_profile.csv`, `middle_vs_outer_summary.csv`, `figures/depth_profile_*.png` |
| 3 | Background fluorescence (explicitly labelled **proxy** unless controls/ROIs exist), with QC flags and a threshold-sensitivity diagnostic | `background_summary.csv`, `threshold_sensitivity.csv`, `figures/background_*.png`, `figures/threshold_sensitivity.png` |
| 4 | Pixel-intensity histograms with identical integer-aligned bins | `intensity_histogram.csv`, `intensity_survival.csv`, `histogram_quantiles.csv`, `figures/histograms_*`, `figures/survival_*` |

`_lng` (deconvolved, full-16-bit-normalised) series are **excluded from every intensity
table** by a code-level guard; they are reserved for morphology QA (not yet implemented).

> **Status:** core pipeline, review gate and 90 synthetic/unit tests are in place.
> **No real microscopy data has been analysed yet.** See [`docs/VALIDATION_STATUS.md`](docs/VALIDATION_STATUS.md)
> for exactly what still needs validation on the actual `.lif`.

## Setup

Requires [uv](https://docs.astral.sh/uv/) (Python ≥ 3.11; tested on 3.13).

```bash
uv sync
```

Dependencies are pinned in `uv.lock`; `uv sync` reproduces the exact environment.
Main libraries: `liffile` (LIF reader), numpy, scipy, scikit-image, pandas, matplotlib, tifffile, PyYAML.

Run the tests (no lab data needed):

```bash
uv run pytest
```

## Usage

### 0. Try the pipeline on synthetic data

```bash
uv run gfap-conc demo
```

Writes `results/synthetic/demo_<timestamp>/`. Synthetic and real outputs are kept apart:

- `demo` only writes under a path containing `synthetic`.
- `review` and `run` refuse such paths.
- Every CSV row carries a `data_origin` column (`synthetic` or `lif`).
- Every synthetic figure is watermarked **SYNTHETIC DATA**.
- `summary.md` and `run_manifest.json` say the data are synthetic.

These numbers are not results.

### 1. Inspect the `.lif` (read-only)

Put the file in `data/raw/` (git-ignored), then:

```bash
uv run gfap-conc inspect data/raw/<file>.lif
```

This writes `results/inspect/inspect_<name>_<timestamp>/`:

- `series_inventory.csv`: one row per series × channel. It covers axes, shape, dtype,
  true bit depth, z-step, XY pixel size, raw min/max, and saturation and zero fractions.
  It also gives a *guessed* raw/`_lng` kind and concentration for each series.
- `acquisition_settings_long.csv`: every captured hardware setting.
- `acquisition_comparability_raw.csv`: settings × raw series, with `differs` flags.
- `lif_xml_header.xml`: the full original metadata, for manual review.
- `input_provenance.json`: path, size and SHA-256 checksum of the file.
- `config_suggestion.yaml`: an **unconfirmed** series list. The GFAP channel is
  deliberately left blank.

### 2. Configure

```bash
cp config/analysis.example.yaml config/analysis.yaml
```

Then edit `config/analysis.yaml`:

- Set `gfap_channel`.
- Paste the series names, concentrations and checksum.
- Leave the analysis defaults unless the inspection shows a reason to change them.

`config/analysis.yaml` is git-ignored because it contains local paths.

### 3. Review tissue mask and z-range (required)

```bash
uv run gfap-conc review --config config/analysis.yaml
```

This writes `results/real/review/review_<timestamp>/` with:

- `figures/review_tissue_mask.png`: the MIPs with the configured mask.
- `figures/review_z_range.png`: an XZ view and per-slice profiles, showing the configured
  range and the auto-detected candidate.
- `review_summary.csv`.
- `review_checklist.md`.

It produces **no** concentration-comparison numbers. Adjust the mask or z-range settings
if needed (e.g. `z_range: manual` with explicit `manual_z_ranges`) and re-run `review`.
When satisfied, record the confirmation in the `review:` block of the config, including
the printed `reviewed_settings_sha256`. That hash covers:

- the input file;
- the channel;
- the series;
- every tissue-mask and z-range setting.

Changing any of them voids the confirmation.

### 4. Run

```bash
uv run gfap-conc run --config config/analysis.yaml
```

`run` refuses to compare concentrations until the review is confirmed for the current
settings. If you need an exploratory look first, pass `--allow-unreviewed`. That run goes
to `preliminary_run_<timestamp>/`, every figure is watermarked **PRELIMINARY**, and the
summary, manifest and QC warnings all say so.

Each run writes a **new** directory: `results/real/run_<timestamp>/`. Existing results are never overwritten. It contains:

- all CSVs listed above, plus `qc_stack_summary.csv` and `qc_warnings.txt`. The QC warnings
  include background-proxy flags (`LARGE_REGION`, `RIGHT_TAIL`, `THRESHOLD_DEPENDENT`,
  `SMALL_REGION`, `EMPTY_REGION`) and a warning if the proxy rises as the antibody is
  diluted;
- `figures/`;
- `summary.md`: auto-generated numbers and QC flags only, ending with a template for the
  human-written interpretation;
- `run_manifest.json`: data origin, review status, git commit, package versions, input
  checksum, resolved config and its hash, thresholds and windows actually used;
- `config_used.yaml`.

The `.lif` is never modified. The run aborts if its checksum does not match the config.

## Project layout

```text
src/gfap_conc/
  model.py        GFAPStack: validated, read-only (Z, Y, X) stack + identity; raw/deconvolved guard
  stats.py        shared descriptive statistics (float64, explicit percentile method)
  projection.py   Analysis 1 - MIP, matched-depth MIP, display ranges
  depth.py        Analysis 2 - per-slice profile, first/middle/last windows
  masks.py        tissue mask, shared foreground threshold (pooled Otsu)
  background.py   Analysis 3 - background proxy / manual ROIs, SBR, CNR
  histograms.py   Analysis 4 - integer-aligned shared bins, survival curves
  pipeline.py     runs all analyses -> tidy tables + QC warnings
  io_lif.py       read-only .lif access (liffile), metadata, inventory
  loading.py      config -> stacks (checksum + channel sanity checks)
  review.py       review figures/checklist, settings hash, confirmation check
  plots.py        figures
  report.py       run manifest + summary.md
  synthetic.py    synthetic stacks for tests/demo
  cli.py          `gfap-conc inspect | review | run | demo`
tests/            unit + end-to-end tests; tests/lif_writer.py builds tiny LIF files
config/           analysis.example.yaml (template)
docs/             methods.md, VALIDATION_STATUS.md, context/ (project background, verbatim)
data/raw/         .lif files (git-ignored)
results/          outputs (git-ignored): inspect/, real/ (review/, run_*, preliminary_run_*), synthetic/
```

## Key rules built into the code

- **Raw intensities are never rescaled.** The per-image p99.8 (`display_vmax_within_image`)
  is a display ceiling only. `normalized_mean_secondary` is labelled as a secondary
  diagnostic because per-image normalisation erases brightness differences by construction.
- **Two display versions of every MIP grid:** one shared scale (an honest brightness
  comparison) and one per-image p99.8 (Rose's style).
- **One foreground threshold, in raw counts, shared by all concentrations.** There are no
  per-image adaptive thresholds.
- **Background is labelled as a proxy** unless it comes from manual ROIs or controls.
  Several QC flags catch a proxy that has probably absorbed faint GFAP, and a
  threshold-sensitivity table shows how background and foreground estimates move with the
  threshold. A rising proxy is never reported as autofluorescence.
- **Background error bars are within-image spread** (robust SD of proxy voxels). They are
  not replicate uncertainty.
- **Tissue mask and z-range must be reviewed on the real images** before concentrations
  are compared.
- **Saturation is measured at the true bit depth** (2^bits − 1, e.g. 4095 for 12-bit data),
  not at 65535.
- **No inferential statistics with n = 1 field per concentration.** Pixels and slices are
  not replicates.

Calculation definitions: [`docs/methods.md`](docs/methods.md).
Step-by-step Windows instructions: [`docs/RUNNING_ON_WINDOWS.md`](docs/RUNNING_ON_WINDOWS.md).
Project background (verbatim copies of the briefing documents): [`docs/context/`](docs/context/).
