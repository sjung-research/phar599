"""Read-only access to Leica .lif files via ``liffile``, plus metadata extraction.

Status: written against the Leica LIF XML layout as parsed by ``liffile`` and tested
only on small synthetic LIF files (tests/lif_writer.py). **Not yet validated on the
real concentration-curve file.** In particular, the hardware-setting attribute names
(detector gain/mode, laser intensity, spectral bands) vary between LAS X versions, so
they are captured generically and must be checked by eye after `gfap-conc inspect`.

Rules enforced here:

* The .lif is only ever opened for reading; nothing is written next to it.
* Axes are taken from metadata labels, never assumed: the GFAP channel is selected
  by index along ``C``, a singleton ``T`` is required, and mosaics / RGB / spectral
  axes raise errors instead of being reshaped.
* Physical sizes are converted to micrometres from each dimension's ``Unit``.
"""

from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass, field
from pathlib import Path
from xml.etree import ElementTree

import numpy as np
import pandas as pd

from .model import GFAPStack, ImageKind

try:  # imported lazily-tolerant so the rest of the package works without it
    from liffile import LifFile
except ImportError:  # pragma: no cover
    LifFile = None  # type: ignore[assignment]

DIM_LABELS = {1: "X", 2: "Y", 3: "Z", 4: "T", 5: "lambda_em", 9: "lambda_ex", 10: "M"}
_UNIT_TO_UM = {"m": 1e6, "mm": 1e3, "um": 1.0, "µm": 1.0, "µm": 1.0, "nm": 1e-3}
DECONVOLVED_TOKEN = "_lng"  # exactly as Rose wrote it; never normalised to another suffix

KEY_SETTINGS = (
    "ObjectiveName", "Magnification", "NumericalAperture", "Zoom", "Pinhole", "PinholeAiry",
    "ScanSpeed", "LineAverage", "Line_Accumulation", "FrameAverage", "FrameAccumulation",
    "Average", "Accumulation", "ScanMode", "IsBidirectional",
)


# ----------------------------------------------------------------- provenance
def sha256_file(path: str | Path, chunk: int = 1 << 22) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while block := fh.read(chunk):
            h.update(block)
    return h.hexdigest()


#: companion files LAS X may write next to a .lif (e.g. extended image data)
COMPANION_SUFFIXES = (".lifext",)


def companion_files(path: str | Path) -> list[Path]:
    """Existing companion files with the same stem (e.g. ``x.lifext`` next to ``x.lif``)."""
    p = Path(path)
    return [p.with_suffix(sfx) for sfx in COMPANION_SUFFIXES if p.suffix.lower() != sfx and p.with_suffix(sfx).exists()]


def file_record(path: str | Path) -> dict:
    p = Path(path)
    return {"path": str(p.resolve()), "sha256": sha256_file(p), "size_bytes": p.stat().st_size}


def _container(image) -> tuple[str, str]:
    parent = getattr(image, "parent", None)
    name = Path(getattr(parent, "filepath", "") or "").name
    kind = getattr(getattr(parent, "type", None), "value", "")
    return name, str(kind)


def _require_liffile() -> None:
    if LifFile is None:  # pragma: no cover
        raise ImportError("liffile is not installed; run `uv sync`")


def _open(path: str | Path):
    _require_liffile()
    return LifFile(str(path), mode="r")


# ------------------------------------------------------------------- metadata
@dataclass
class DimensionMeta:
    label: str
    dim_id: int
    n: int
    length_um: float | None  # signed, first to last element
    unit: str

    @property
    def step_um(self) -> float | None:
        if self.length_um is None or self.n < 2:
            return None
        return abs(self.length_um) / (self.n - 1)


@dataclass
class ChannelMeta:
    index: int  # position along liffile's C axis (sorted by BytesInc)
    name: str
    resolution_bits: int
    data_type: int
    lut_name: str


@dataclass
class SeriesMeta:
    path: str
    name: str
    dims: tuple[str, ...]
    sizes: dict[str, int]
    dtype: str
    dimensions: dict[str, DimensionMeta]
    channels: list[ChannelMeta]
    hardware: dict[str, str] = field(default_factory=dict)
    first_timestamp: str | None = None
    has_tilescan: bool = False

    @property
    def z_step_um(self) -> float | None:
        d = self.dimensions.get("Z")
        return d.step_um if d else None

    @property
    def pixel_size_um(self) -> tuple[float, float] | None:
        y, x = self.dimensions.get("Y"), self.dimensions.get("X")
        if y is None or x is None or y.step_um is None or x.step_um is None:
            return None
        return (y.step_um, x.step_um)


def _to_um(value: float, unit: str) -> float | None:
    f = _UNIT_TO_UM.get(unit.strip())
    return value * f if f is not None else None


def parse_dimensions(xml: ElementTree.Element) -> dict[str, DimensionMeta]:
    out = {}
    for d in xml.findall("./Data/Image/ImageDescription/Dimensions/DimensionDescription"):
        dim_id = int(d.attrib["DimID"])
        label = DIM_LABELS.get(dim_id, f"dim{dim_id}")
        unit = d.attrib.get("Unit", "")
        out[label] = DimensionMeta(label, dim_id, int(d.attrib["NumberOfElements"]),
                                   _to_um(float(d.attrib.get("Length", "nan")), unit), unit)
    return out


def parse_channels(xml: ElementTree.Element, names: list[str] | None = None) -> list[ChannelMeta]:
    els = sorted(
        xml.findall("./Data/Image/ImageDescription/Channels/ChannelDescription"),
        key=lambda e: int(e.attrib.get("BytesInc", 0)),
    )
    out = []
    for i, e in enumerate(els):
        out.append(ChannelMeta(i, names[i] if names and i < len(names) else f"Ch{i}",
                               int(e.attrib["Resolution"]), int(e.attrib.get("DataType", 0)),
                               e.attrib.get("LUTName", "")))
    return out


def flatten_hardware(xml: ElementTree.Element) -> dict[str, str]:
    """All attributes of the acquisition-relevant HardwareSetting elements, flattened.

    Keys look like ``ATLConfocalSettingDefinition[0].Zoom`` or
    ``Detector[HyD S 3].Gain``. Inactive detectors are skipped.
    """
    hs = xml.find('./Data/Image/Attachment[@Name="HardwareSetting"]')
    out: dict[str, str] = {}
    if hs is None:
        return out
    counters: dict[str, int] = {}
    for el in hs.iter():
        tag = el.tag
        if tag == "ATLConfocalSettingDefinition":
            i = counters.get(tag, 0)
            counters[tag] = i + 1
            for k, v in el.attrib.items():
                out[f"{tag}[{i}].{k}"] = v
        elif tag == "Detector":
            if el.attrib.get("IsActive", "1") not in ("1", "true", "True"):
                continue
            key = el.attrib.get("Name", f"#{counters.setdefault('Detector', 0)}")
            for k, v in el.attrib.items():
                out[f"Detector[{key}].{k}"] = v
        elif tag in ("LaserLineSetting", "MultiBand"):
            ident = el.attrib.get("LaserLine") or el.attrib.get("Channel") or str(counters.get(tag, 0))
            counters[tag] = counters.get(tag, 0) + 1
            for k, v in el.attrib.items():
                out[f"{tag}[{ident}].{k}"] = v
    return out


def key_settings(hardware: dict[str, str]) -> dict[str, str]:
    """Selected scan settings of the first confocal setting block, if present."""
    return {k: hardware[f"ATLConfocalSettingDefinition[0].{k}"]
            for k in KEY_SETTINGS if f"ATLConfocalSettingDefinition[0].{k}" in hardware}


def describe_image(image) -> SeriesMeta:
    xml = image.xml_element
    names = None
    coords = getattr(image, "coords", {}) or {}
    if "C" in coords:
        names = [str(c) for c in coords["C"]]
    ts = image.timestamps if hasattr(image, "timestamps") else None
    return SeriesMeta(
        path=image.path,
        name=image.name,
        dims=tuple(image.dims),
        sizes=dict(image.sizes),
        dtype=str(image.dtype),
        dimensions=parse_dimensions(xml),
        channels=parse_channels(xml, names),
        hardware=flatten_hardware(xml),
        first_timestamp=str(ts[0]) if ts is not None and len(ts) else None,
        has_tilescan=getattr(image, "tilescan", None) is not None,
    )


def list_series(path: str | Path) -> list[SeriesMeta]:
    with _open(path) as lif:
        return [describe_image(im) for im in lif.images]


# ------------------------------------------------------------ name guessing
_CONC_RES = (
    re.compile(r"(?<![\d.])1\s*[:_\-]\s*(250|500|1000|2000)(?!\d)"),
    re.compile(r"(?<![\d.])(250|500|1000|2000)(?!\d)"),
)


def guess_concentration(name: str) -> str | None:
    """Proposed concentration from a series name. A guess only - confirm in config."""
    for rx in _CONC_RES:
        hits = {m.group(1) for m in rx.finditer(name)}
        if len(hits) == 1:
            return f"1:{hits.pop()}"
        if len(hits) > 1:
            return None
    return None


def guess_kind(name: str) -> ImageKind:
    return ImageKind.DECONVOLVED if DECONVOLVED_TOKEN in name else ImageKind.RAW


# ----------------------------------------------------------------- pixel I/O
_ALLOWED_DIMS = {"Z", "Y", "X", "C"}


def channel_zyx(array: np.ndarray, dims: tuple[str, ...], channel: int, label: str = "") -> np.ndarray:
    """Select one channel and return an explicit (Z, Y, X) array.

    ``dims`` are liffile's labels for ``array`` (singleton axes already squeezed).
    """
    dims = tuple(dims)
    if array.ndim != len(dims):
        raise ValueError(f"{label}: array ndim {array.ndim} does not match dims {dims}")
    extra = set(dims) - _ALLOWED_DIMS
    if extra:
        raise ValueError(
            f"{label}: unsupported axes {sorted(extra)} in {dims} (time series, mosaic, RGB or spectral data); "
            "handle explicitly before quantification"
        )
    if "C" in dims:
        n_c = array.shape[dims.index("C")]
        if not 0 <= channel < n_c:
            raise ValueError(f"{label}: channel {channel} out of range (C={n_c})")
        array = np.take(array, channel, axis=dims.index("C"))
        dims = tuple(d for d in dims if d != "C")
    elif channel != 0:
        raise ValueError(f"{label}: single-channel image, but channel {channel} requested")
    if "Z" not in dims:
        raise ValueError(f"{label}: no Z axis - not a z-stack (dims {dims})")
    order = [dims.index(d) for d in ("Z", "Y", "X")]
    return np.ascontiguousarray(np.transpose(array, order))


def read_channel(path: str | Path, series_path: str, channel: int) -> tuple[np.ndarray, SeriesMeta]:
    with _open(path) as lif:
        image = lif.images[series_path]
        meta = describe_image(image)
        arr = channel_zyx(image.asarray(), meta.dims, channel, series_path)
    return arr, meta


def load_stack(
    path: str | Path,
    series_path: str,
    channel: int,
    concentration: str,
    kind: ImageKind | str,
    bit_depth_override: int | None = None,
) -> GFAPStack:
    arr, meta = read_channel(path, series_path, channel)
    bits = bit_depth_override or meta.channels[channel].resolution_bits
    return GFAPStack(
        series_name=series_path,
        concentration=concentration,
        kind=ImageKind(kind),
        channel_index=channel,
        data=arr,
        bit_depth=bits,
        z_step_um=meta.z_step_um,
        pixel_size_um=meta.pixel_size_um,
        intensity_unit="raw counts",
        source=f"{Path(path).name}::{series_path}",
    )


# ------------------------------------------------------------------ inventory
def inventory(path: str | Path, pixel_stats: bool = True) -> tuple[pd.DataFrame, pd.DataFrame]:
    """One row per (series, channel), and a long table of all hardware settings."""
    rows, settings = [], []
    with _open(path) as lif:
        for image in lif.images:
            meta = describe_image(image)
            arr = image.asarray() if pixel_stats else None
            z = meta.dimensions.get("Z")
            for ch in meta.channels:
                container_file, container_type = _container(image)
                row = {
                    "series_path": meta.path,
                    "series_name": meta.name,
                    "container_file": container_file,
                    "container_type": container_type,
                    "image_kind_guess": guess_kind(meta.path).value,
                    "concentration_guess": guess_concentration(meta.name),
                    "channel_index": ch.index,
                    "channel_name": ch.name,
                    "channel_lut": ch.lut_name,
                    "dims": "".join(meta.dims),
                    "shape": "x".join(str(meta.sizes[d]) for d in meta.dims),
                    "Z": meta.sizes.get("Z", 1),
                    "Y": meta.sizes.get("Y", 1),
                    "X": meta.sizes.get("X", 1),
                    "T": meta.sizes.get("T", 1),
                    "C": meta.sizes.get("C", 1),
                    "dtype": meta.dtype,
                    "bit_depth": ch.resolution_bits,
                    "z_step_um": meta.z_step_um,
                    "z_extent_um": abs(z.length_um) if z and z.length_um is not None else None,
                    "z_direction": (None if not z or not z.length_um else ("+" if z.length_um > 0 else "-")),
                    "pixel_size_y_um": meta.pixel_size_um[0] if meta.pixel_size_um else None,
                    "pixel_size_x_um": meta.pixel_size_um[1] if meta.pixel_size_um else None,
                    "first_timestamp": meta.first_timestamp,
                    "tilescan": meta.has_tilescan,
                    **{f"setting_{k}": v for k, v in key_settings(meta.hardware).items()},
                }
                if arr is not None:
                    try:
                        data = channel_zyx(arr, meta.dims, ch.index, meta.path)
                        sat = 2**ch.resolution_bits - 1
                        n = data.size
                        row.update(
                            raw_min=int(data.min()), raw_max=int(data.max()), saturation_value=sat,
                            saturation_fraction=float(np.count_nonzero(data >= sat)) / n,
                            zero_fraction=float(np.count_nonzero(data == 0)) / n,
                            values_exceed_bit_depth=bool(data.max() > sat),
                        )
                    except ValueError as exc:
                        row["pixel_stats_error"] = str(exc)
                rows.append(row)
            settings += [{"series_path": meta.path, "key": k, "value": v} for k, v in meta.hardware.items()]
    return pd.DataFrame(rows), pd.DataFrame(settings, columns=["series_path", "key", "value"])


def comparability_table(settings_long: pd.DataFrame, series_paths: list[str]) -> pd.DataFrame:
    """Settings as rows x series as columns, with a ``differs`` flag."""
    sub = settings_long[settings_long["series_path"].isin(series_paths)]
    if sub.empty:
        return pd.DataFrame(columns=["key", "differs"])
    wide = sub.pivot_table(index="key", columns="series_path", values="value", aggfunc="first")
    wide = wide.reindex(columns=[p for p in series_paths if p in wide.columns])
    wide["differs"] = wide.apply(lambda r: r.dropna().nunique() > 1 or r.isna().any(), axis=1)
    return wide.reset_index().sort_values(["differs", "key"], ascending=[False, True])


def xml_header(path: str | Path) -> str:
    with _open(path) as lif:
        return lif.xml_header()
