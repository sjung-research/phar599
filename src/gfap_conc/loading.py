"""Turn a validated :class:`ProjectConfig` into GFAP stacks read from the .lif."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from .config import ProjectConfig
from .io_lif import companion_files, file_record, list_series, load_stack, read_channel, sha256_file
from .model import GFAPStack, ImageKind
from .projection import max_projection

_NUCLEAR_RE = re.compile(r"dapi|hoechst|405", re.IGNORECASE)


class ConfigIncompleteError(ValueError):
    pass


def resolve_lif_path(cfg: ProjectConfig, config_dir: Path) -> Path:
    if not cfg.input.lif_path:
        raise ConfigIncompleteError("input.lif_path is not set")
    p = Path(cfg.input.lif_path)
    return p if p.is_absolute() else (config_dir / p).resolve()


def verify_checksum(path: Path, expected: str | None) -> str:
    actual = sha256_file(path)
    if expected and actual.lower() != expected.lower():
        raise ValueError(f"SHA-256 mismatch for {path.name}: expected {expected}, got {actual}")
    return actual


def load_configured_stacks(
    cfg: ProjectConfig, config_dir: Path
) -> tuple[list[GFAPStack], dict[str, np.ndarray], list[str], dict]:
    """Return (stacks, reference MIPs for tissue masking, warnings, provenance)."""
    if cfg.input.gfap_channel is None:
        raise ConfigIncompleteError("input.gfap_channel must be set after inspecting the file")
    if not cfg.series:
        raise ConfigIncompleteError("no series configured; copy them from the inspect output")
    lif = resolve_lif_path(cfg, config_dir)
    digest = verify_checksum(lif, cfg.input.lif_sha256)
    warnings = []
    if not cfg.input.lif_sha256:
        warnings.append("input.lif_sha256 not set in config; input file identity was not verified")

    stacks, refs = [], {}
    skipped = [s.name for s in cfg.series if s.kind is not ImageKind.RAW]
    if skipped:
        # Deconvolved images never enter intensity analyses, and morphology QA is not
        # implemented yet, so their pixels are not loaded (saves memory).
        warnings.append("Configured deconvolved series not loaded and excluded from all intensity analyses: "
                        + ", ".join(skipped))
    for spec in cfg.series:
        if spec.kind is not ImageKind.RAW:
            continue
        st = load_stack(lif, spec.name, cfg.input.gfap_channel, spec.concentration, spec.kind,
                        cfg.input.bit_depth_override, cfg.input.intensity_unit)
        stacks.append(st)
        if cfg.analysis.tissue_mask.reference == "reference_channel_mip":
            if cfg.input.reference_channel is None:
                raise ConfigIncompleteError("tissue mask uses reference_channel_mip but input.reference_channel is unset")
            ref, _ = read_channel(lif, spec.name, cfg.input.reference_channel)
            refs[spec.name] = max_projection(ref)

    # Guard against quantifying the nuclear channel as GFAP.
    metas = {m.path: m for m in list_series(lif)}
    for spec in cfg.series:
        m = metas.get(spec.name)
        if m is None:
            continue
        ch = m.channels[cfg.input.gfap_channel]
        evidence = f"{ch.name} {ch.lut_name}"
        if _NUCLEAR_RE.search(ch.name):
            raise ValueError(f"{spec.name}: channel {ch.index} is named {ch.name!r}; this looks like DAPI, not GFAP")
        if ch.lut_name.lower() == "blue":
            warnings.append(f"{spec.name}: GFAP channel {ch.index} uses a blue LUT ({evidence}); confirm it is AF647")
    provenance = {"lif_path": str(lif), "lif_sha256": digest, "lif_size_bytes": lif.stat().st_size,
                  "companion_files": [file_record(c) for c in companion_files(lif)]}
    return stacks, refs, warnings, provenance
