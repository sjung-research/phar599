# Session log: building and running the GFAP concentration-curve pipeline

Covers one Claude Code cloud session, 2026-11-20 to 2026-11-23, on branch
`claude/blissful-cori-wjn019` of `sjung-research/phar599`.

## Authorship note (read first)

The Python pipeline in this repository was **written by Claude (an AI coding assistant) in
this session**, at Sewon's direction and with Sewon reviewing each step. Sewon supplied:
- the project context (three briefing documents, kept verbatim in `docs/context/`);
- the requirements and their priorities;
- all decisions the tool could not make itself: channel, series mapping, and confirmation
  of the tissue mask and z-range;
- every run on the real data.

Rose suggested that Sewon first build the analysis independently and then compare it with
her earlier script. This code should therefore be presented to her as **AI-assisted**, not
as an unaided implementation. The design choices and their reasons are recorded below and
in `docs/methods.md`, so they can be compared with hers point by point.

## Where things ran and what Claude saw

- The raw `.lif` (2.40 GB) and `.lifext` (0.38 GB) **never left Sewon's laptop**. Claude
  wrote and tested code in a cloud container, Sewon pulled it with git and ran it locally
  (Windows 11, PowerShell, uv).
- Claude saw only what Sewon pasted into the chat:
  - console output;
  - the series inventory and acquisition settings;
  - the review checklist and summary CSV;
  - `summary.md`, `qc_warnings.txt`;
  - the figures (PNG).
- Run outputs (`results/`) are git-ignored and exist only on the laptop.

## Timeline

| Date | Phase | Commits |
|---|---|---|
| 11-20 | Read the three briefing documents; proposed tech stack, workflow, pitfalls and implementation order (no code). | – |
| 11-20 | Built the core package, Rose's four analyses, synthetic tests and docs. | `aefc639` |
| 11-20 | Review requests from Sewon: background QC flags, threshold sensitivity, mandatory review gate, voxel/pixel wording, label overlap, error-bar caption, synthetic vs real separation. | `1b06834` |
| 11-22 | Laptop setup. Found a `.lifext` companion file; added provenance for it and a Windows guide. | `778bde9` |
| 11-22 | Two tests failed on Windows (Korean cp949 locale); all file I/O made explicit UTF-8. | `f4ebd9d` |
| 11-22 | First `inspect` of the real file: the suffix is `_Lng` (not `_lng`), and raw maxima are nearly identical. | `7e96faf` |
| 11-22 | Acquisition settings checked: identical across the four images, photon counting. Unit label added; excluded `_Lng` pixels no longer loaded. | `fdd5bb5` |
| 11-22 | First review on real images: the automatic z-range rule failed and was fixed; XZ view made readable. | `8666cb5` |
| 11-22 | Second review: tissue mask `none`, z-range `full` confirmed by Sewon. | `e7ccd64` |
| 11-22 | First confirmed run. Histogram bins too coarse and `LARGE_REGION` flag miscalibrated; both fixed. | `e5552bf` |
| 11-22 | Re-run with fine bins revealed soft detector compression; recorded. | `811e1b7` |
| 11-27 | Compared outputs with Rose's first-staining PDFs (and a ChatGPT review of them). Added figures in Rose's layout: MIP with its pixel histogram and p99.8 line; depth slice gallery at 0/25/50/75/90% using her slice rule. | `b491535` |
| 11-23 | English summary report for Rose (Claude Docs, docx and md in `reports/`), this log, and `docs/AI_HANDOFF.md`. | latest commits |

## What was built

A Python package, `gfap-conc` (`src/gfap_conc/`), managed with uv (`pyproject.toml`,
`uv.lock`, Python 3.13).

| Module | Role |
|---|---|
| `model.py` | `GFAPStack`: validated, read-only `(Z, Y, X)` array with concentration, series name, channel, bit depth, units; raw vs deconvolved guard |
| `stats.py` | Shared statistics (float64, explicit percentile method) |
| `projection.py` | Analysis 1: MIP, matched-depth MIP, display ranges |
| `depth.py` | Analysis 2: per-slice profile, first/middle/last windows, z-range rule |
| `masks.py` | Tissue mask, one shared foreground threshold (pooled Otsu) |
| `background.py` | Analysis 3: background proxy, ROIs, SBR/CNR, QC flags, threshold sensitivity |
| `histograms.py` | Analysis 4: integer-aligned shared bins, survival curves |
| `pipeline.py` | Runs all analyses; tidy tables; QC warnings |
| `io_lif.py`, `loading.py` | Read-only `.lif` access (liffile), metadata, inventory, checksum |
| `review.py` | Review figures and checklist; settings hash; confirmation check |
| `plots.py`, `report.py` | Figures; `summary.md` and `run_manifest.json` |
| `synthetic.py` | Synthetic stacks for tests and the demo |
| `cli.py` | `gfap-conc inspect`, `review`, `run`, `demo` |

There are 90 tests (`tests/`). They use synthetic stacks and small LIF files written to the
documented format by `tests/lif_writer.py`; no lab data is required.

## Key design decisions

- **Raw images only.** `_Lng` images (deconvolved, rescaled to 16 bit) are blocked from every
  intensity function, as Rose specified.
- **Raw means and p99.8 are kept separate.** p99.8 is a display ceiling only. Per-image
  normalisation is reported as a secondary column, because it erases brightness differences
  between images by construction.
- **One shared foreground threshold** for all dilutions, plus a sensitivity sweep, instead
  of per-image thresholds.
- **Background is always labelled a proxy** (no no-primary controls), with QC flags, and is
  never reported as autofluorescence.
- **A mandatory human review** of tissue mask and z-range is required before any comparison.
  It is tied to a hash of the reviewed settings.
- **Provenance:** the input SHA-256 is checked on every run; outputs never overwrite earlier
  ones; each run writes a manifest with code commit, package versions and the resolved
  config.
- **Descriptive statistics only.** There is one field per concentration, so no p-values.

## Problems found and fixed during the session

| Problem | How it surfaced | Fix |
|---|---|---|
| Histogram bins could fall entirely above the maximum value | Synthetic test | Emit only the bins needed to cover the range |
| A shared threshold biases the background proxy upward for weak stains | Synthetic test | Kept the design; added QC flags, sensitivity sweep, documentation |
| Raw depth ratios include the background floor, which pulls them toward 1 | Synthetic demo | Documented; background-subtracted variant deferred |
| Tests failed on Korean Windows (cp949) | Sewon's laptop | Explicit UTF-8 everywhere; BOM-tolerant config loading |
| Sewon pasted YAML into PowerShell | Laptop | Config written by a single PowerShell command instead |
| Deconvolved suffix is `_Lng`, so it was guessed as raw | Real `inspect` | Case-insensitive match with exact spelling reported |
| `.lifext` images not read through the `.lif` | Real `inspect` | Recorded with checksum; not needed for the analysis |
| Automatic z-range rule cut the stacks at their mid-depth dip | First real review | Level measured from background floor; only end slices trimmed |
| XZ review view was a thin strip | First real review | z stretched, factor shown |
| Histogram bins 256 counts wide, data ≤ 1,429 | First real run | Bins span the observed range (6 counts) |
| `LARGE_REGION` flagged every stack (GFAP fills only 2–4% of voxels) | First real run | Limit 0.85 → 0.97; no result value affected |

## Real-data results in brief

Full report: *Mouse anti-GFAP concentration curve: first quantitative report*
(<https://claude.ai/artifact/EMni7yBiV8BouX5G2bJP55>). The link is private until shared from
the page's Share menu.

- Acquisition settings were identical:
  - 20x/0.75 dry objective, zoom 2, 0.132 µm pixels;
  - HyD S detectors in photon counting, gain 2.5;
  - 638 nm laser at 2.5%;
  - line accumulation 4, frame average 2.
- The four images come from four different sections.
- **No dilution effect** on MIP mean (219.5, 204.6, 288.3, 199.4 counts), depth retention,
  background (21–28 counts) or distribution shape. 1:1000 is higher because its field
  contains more GFAP-positive tissue, not because its structures are brighter.
- **The mid-depth dip is the same at every dilution:** about 55–63% of the near-surface peak
  for the mean.
- **Soft detector saturation near 1,400 counts:** maxima and p99.8 measure the ceiling, so
  the comparison relies on means, medians and the 100–1,200-count range.

## Open items

- Linearity test at lower excitation or accumulation; then re-image below the ceiling.
- Several fields per dilution in matched regions; extend to 1:4000 and 1:8000; no-primary
  controls.
- Confirm region comparability of the four sections with Rose.
- Not yet implemented: `_Lng` morphology QA, background-subtracted depth ratios, handling of
  several fields per concentration, Bio-Formats cross-check.
- Compare this pipeline with Rose's earlier script, which has not been seen yet.

## How to reproduce

See `README.md` and `docs/RUNNING_ON_WINDOWS.md`. In short:

```bash
uv sync
uv run pytest
uv run gfap-conc inspect "<path>\Sewon_GFAP_curve.lif"
uv run gfap-conc review --config config\analysis.yaml
uv run gfap-conc run --config config\analysis.yaml
```

`config/analysis.yaml` is git-ignored. Its content for this run is recorded in each run's
`config_used.yaml` and `run_manifest.json` on the laptop.
