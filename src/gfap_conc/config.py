"""YAML configuration -> typed, validated settings.

Unknown keys are rejected so typos (e.g. ``exclusion_radius`` vs
``exclusion_radius_px``) fail loudly instead of silently using a default.
"""

from __future__ import annotations

import dataclasses
import hashlib
import typing
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

from .model import ImageKind, parse_dilution


@dataclass
class TissueMaskSettings:
    method: str = "none"  # none | absolute | otsu
    reference: str = "gfap_mip"  # gfap_mip | reference_channel_mip
    threshold: float | None = None
    sigma_px: float = 2.0
    min_object_px: int = 500


@dataclass
class ForegroundSettings:
    #: number (raw counts), "pooled_otsu", or null (skip foreground-based metrics)
    threshold: float | str | None = "pooled_otsu"
    sample_per_stack: int = 200_000
    seed: int = 0
    #: diagnostic only: multiples of the resolved threshold at which background/foreground are re-estimated
    sensitivity_factors: list[float] = field(default_factory=lambda: [0.5, 0.75, 1.0, 1.25, 1.5, 2.0])


@dataclass
class BackgroundSettings:
    exclusion_radius_px: int = 5
    # --- QC limits for the background proxy (warnings only; nothing is auto-corrected) ---
    #: warn if the proxy covers more than this fraction of tissue voxels
    qc_max_fraction_of_tissue: float = 0.85
    #: warn if (p95 - median) / (median - p05) of proxy voxels exceeds this (right tail = faint signal)
    qc_max_tail_asymmetry: float = 2.0
    #: warn if the proxy median changes by more than this fraction between the sensitivity
    #: factors bracketing 1.0 (default 0.75x and 1.25x of the threshold)
    qc_max_threshold_dependence: float = 0.10
    #: half-open (y0, y1, x0, x1) rectangles in original pixel coordinates, per series name
    manual_rois: dict[str, list[list[int]]] = field(default_factory=dict)


@dataclass
class DepthSettings:
    z_range: str = "full"  # full | auto | manual
    #: for z_range: manual -> {series name: [first, last]} inclusive slice indices, chosen during review
    manual_z_ranges: dict[str, list[int]] = field(default_factory=dict)
    auto_metric: str = "p99_8"
    auto_rel_level: float = 0.5
    window_fraction: float = 0.2
    retention_metrics: list[str] = field(default_factory=lambda: ["mean", "median", "p99_8"])


@dataclass
class HistogramSettings:
    n_bins: int = 256
    zoom_max: float | None = None  # null -> 1.2 x largest per-image p99.8 of the MIPs


@dataclass
class AnalysisSettings:
    display_percentile: float = 99.8
    tissue_mask: TissueMaskSettings = field(default_factory=TissueMaskSettings)
    foreground: ForegroundSettings = field(default_factory=ForegroundSettings)
    background: BackgroundSettings = field(default_factory=BackgroundSettings)
    depth: DepthSettings = field(default_factory=DepthSettings)
    histogram: HistogramSettings = field(default_factory=HistogramSettings)


@dataclass
class SeriesSpec:
    name: str  # exact image path inside the .lif, as listed by `gfap-conc inspect`
    concentration: str  # "1:250" ...
    kind: ImageKind = ImageKind.RAW
    pair_of: str | None = None  # for deconvolved: name of the matching raw series

    def __post_init__(self) -> None:
        self.kind = ImageKind(self.kind)
        parse_dilution(self.concentration)


@dataclass
class InputSettings:
    lif_path: str = ""
    lif_sha256: str | None = None  # verified before every run when set
    gfap_channel: int | None = None  # REQUIRED; index in liffile's C order
    reference_channel: int | None = None  # e.g. DAPI, only for tissue masking
    bit_depth_override: int | None = None


@dataclass
class ReviewSettings:
    """Human confirmation that tissue mask and z-range were checked on the real images
    (`gfap-conc review`). `run` refuses to compare concentrations until both are confirmed
    for exactly the current settings (matched by hash)."""

    tissue_mask_confirmed: bool = False
    z_range_confirmed: bool = False
    reviewed_settings_sha256: str | None = None
    reviewed_by: str | None = None
    notes: str = ""


@dataclass
class ProjectConfig:
    input: InputSettings = field(default_factory=InputSettings)
    series: list[SeriesSpec] = field(default_factory=list)
    analysis: AnalysisSettings = field(default_factory=AnalysisSettings)
    review: ReviewSettings = field(default_factory=ReviewSettings)
    output_dir: str = "results/real"


def _build(cls: type, data: Any, where: str) -> Any:
    if not dataclasses.is_dataclass(cls):
        return data
    if data is None:
        return cls()
    if not isinstance(data, dict):
        raise ValueError(f"{where}: expected a mapping, got {type(data).__name__}")
    hints = typing.get_type_hints(cls)
    names = {f.name for f in dataclasses.fields(cls)}
    unknown = set(data) - names
    if unknown:
        raise ValueError(f"{where}: unknown key(s) {sorted(unknown)}; allowed: {sorted(names)}")
    kwargs = {}
    for key, value in data.items():
        hint = hints[key]
        origin = typing.get_origin(hint)
        args = typing.get_args(hint)
        if dataclasses.is_dataclass(hint):
            kwargs[key] = _build(hint, value, f"{where}.{key}")
        elif origin is list and args and dataclasses.is_dataclass(args[0]):
            kwargs[key] = [_build(args[0], v, f"{where}.{key}[{i}]") for i, v in enumerate(value or [])]
        else:
            kwargs[key] = value
    return cls(**kwargs)


def load_config(path: str | Path) -> ProjectConfig:
    path = Path(path)
    with path.open(encoding="utf-8") as fh:
        raw = yaml.safe_load(fh) or {}
    cfg: ProjectConfig = _build(ProjectConfig, raw, path.name)
    validate_config(cfg)
    return cfg


def validate_config(cfg: ProjectConfig) -> None:
    a = cfg.analysis
    if a.tissue_mask.method not in ("none", "absolute", "otsu"):
        raise ValueError(f"analysis.tissue_mask.method: invalid {a.tissue_mask.method!r}")
    if a.tissue_mask.reference not in ("gfap_mip", "reference_channel_mip"):
        raise ValueError(f"analysis.tissue_mask.reference: invalid {a.tissue_mask.reference!r}")
    if a.depth.z_range not in ("full", "auto", "manual"):
        raise ValueError(f"analysis.depth.z_range: invalid {a.depth.z_range!r}")
    for name, zr in a.depth.manual_z_ranges.items():
        if len(zr) != 2 or not 0 <= int(zr[0]) <= int(zr[1]):
            raise ValueError(f"analysis.depth.manual_z_ranges[{name!r}]: need [first, last] with 0 <= first <= last")
    if a.depth.z_range == "manual":
        missing = [s.name for s in cfg.series if s.kind is ImageKind.RAW and s.name not in a.depth.manual_z_ranges]
        if missing:
            raise ValueError(f"analysis.depth.z_range is 'manual' but no range is given for {missing}")
    f = a.foreground.sensitivity_factors
    if not f or any(v <= 0 for v in f) or 1.0 not in f:
        raise ValueError("analysis.foreground.sensitivity_factors must be positive and include 1.0")
    thr = a.foreground.threshold
    if isinstance(thr, str) and thr != "pooled_otsu":
        raise ValueError(f"analysis.foreground.threshold: invalid {thr!r}")
    if not 0 < a.display_percentile <= 100:
        raise ValueError("analysis.display_percentile must be in (0, 100]")
    names = [s.name for s in cfg.series]
    if len(names) != len(set(names)):
        raise ValueError("series: duplicate series names")
    raw_names = {s.name for s in cfg.series if s.kind is ImageKind.RAW}
    for s in cfg.series:
        if s.kind is ImageKind.DECONVOLVED and s.pair_of is not None and s.pair_of not in raw_names:
            raise ValueError(f"series {s.name!r}: pair_of {s.pair_of!r} is not a configured raw series")


def config_sha256(cfg: ProjectConfig) -> str:
    """Hash of the *resolved* configuration (defaults included)."""
    blob = yaml.safe_dump(_to_plain(dataclasses.asdict(cfg)), sort_keys=True).encode()
    return hashlib.sha256(blob).hexdigest()


def _to_plain(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {k: _to_plain(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_to_plain(v) for v in obj]
    if isinstance(obj, ImageKind):
        return obj.value
    return obj


def resolved_config_dict(cfg: ProjectConfig) -> dict:
    return _to_plain(dataclasses.asdict(cfg))
