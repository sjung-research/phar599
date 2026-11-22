# Running on a Windows laptop (PowerShell)

The `.lif` stays where it is (e.g. on an external SSD). It is never copied into the
repository and never modified.

## One-time setup

1. Install **Git for Windows** (https://git-scm.com/download/win). Keep the defaults; Git
   Credential Manager will open a browser to log in to GitHub on the first clone.
2. Install **uv** in PowerShell:

   ```powershell
   powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
   ```

   Close and reopen PowerShell so `uv` is on the PATH. `uv` downloads the right Python
   (3.13) by itself.
3. Clone the analysis branch and install the locked environment:

   ```powershell
   cd $HOME\Documents
   git clone -b claude/blissful-cori-wjn019 https://github.com/sjung-research/phar599.git
   cd phar599
   uv sync
   uv run pytest
   ```

   All tests should pass before touching real data.

## Protect the raw file (recommended)

```powershell
attrib +R "E:\2026NovemberGFAPlif\Sewon_GFAP_curve.lif"
attrib +R "E:\2026NovemberGFAPlif\Sewon_GFAP_curve.lifext"
```

This sets the read-only attribute. The pipeline only opens files for reading anyway; this
also guards against other software (LAS X, Fiji) saving over them.

## Step 1: inspect (read-only)

```powershell
uv run gfap-conc inspect "E:\2026NovemberGFAPlif\Sewon_GFAP_curve.lif"
```

Output goes to `results\inspect\inspect_Sewon_GFAP_curve_lif_<timestamp>\`. If the
command prints a NOTE about the `.lifext` file, inspect that file too:

```powershell
uv run gfap-conc inspect "E:\2026NovemberGFAPlif\Sewon_GFAP_curve.lifext"
```

If it runs out of memory or is very slow, add `--no-pixel-stats` (metadata only).

**Stop here** and check `series_inventory.csv`, `acquisition_comparability_raw.csv` and
`config_suggestion.yaml` before configuring anything.

## Step 2: configure

```powershell
Copy-Item config\analysis.example.yaml config\analysis.yaml
notepad config\analysis.yaml
```

In YAML, write Windows paths in **single quotes** so backslashes are kept literally:

```yaml
input:
  lif_path: 'E:\2026NovemberGFAPlif\Sewon_GFAP_curve.lif'
  lif_sha256: <from config_suggestion.yaml>
  gfap_channel: <confirmed AF647 channel index>
```

Copy the `series` entries from `config_suggestion.yaml`. Fix every concentration and kind,
and delete series that are not part of the curve. `config/analysis.yaml` is git-ignored.

## Step 3: review tissue mask and z-range

```powershell
uv run gfap-conc review --config config\analysis.yaml
```

Open `review_checklist.md` and the two figures in `results\real\review\review_<timestamp>\`.
When satisfied, fill in the `review:` block of the config with the printed hash.

## Step 4: run

```powershell
uv run gfap-conc run --config config\analysis.yaml
```

Results: `results\real\run_<timestamp>\` (CSVs, `figures\`, `summary.md`,
`run_manifest.json`).

## Getting updates

```powershell
git pull
uv sync
```
