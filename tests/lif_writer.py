"""Minimal writer for LIF version-2 files, used ONLY to test the reader adapter.

Layout (Leica LIF spec v3.2, as parsed by liffile):
  header  : <IIBI> 0x70, 2*nchars+5, 0x2A, nchars, then UTF-16LE XML
  blocks  : <IIBQBI> 0x70, 0, 0x2A, nbytes, 0x2A, nchars, UTF-16LE block id, raw bytes

This mimics the structure LAS X writes but is not LAS X output: hardware-setting
attribute names here are illustrative. Passing tests prove our code agrees with
liffile on this structure - not that it handles every real LAS X file.
"""

from __future__ import annotations

import struct
from dataclasses import dataclass, field
from pathlib import Path
from xml.sax.saxutils import quoteattr

import numpy as np


@dataclass
class LifTestImage:
    name: str
    data: np.ndarray  # (Z, C, Y, X) uint16
    resolution_bits: int = 12
    voxel_um: tuple[float, float, float] = (0.5, 0.2, 0.2)  # z, y, x
    luts: tuple[str, ...] = ("Red", "Blue")
    hardware_attrs: dict[str, str] = field(default_factory=dict)
    detectors: list[dict[str, str]] = field(default_factory=list)


def _image_xml(img: LifTestImage, block_id: str) -> str:
    nz, nc, ny, nx = img.data.shape
    item = 2
    chans = "".join(
        f'<ChannelDescription DataType="0" ChannelTag="0" Resolution="{img.resolution_bits}" '
        f'NameOfMeasuredQuantity="" Min="0" Max="{2**img.resolution_bits - 1}" Unit="" '
        f'LUTName="{img.luts[c] if c < len(img.luts) else "Gray"}" IsLUTInverted="0" '
        f'BytesInc="{c * nx * ny * item}" BitInc="0"/>'
        for c in range(nc)
    )
    dz, dy, dx = (v * 1e-6 for v in img.voxel_um)  # LAS X stores metres

    def dim(i, n, step, inc):
        return (f'<DimensionDescription DimID="{i}" NumberOfElements="{n}" Origin="0" '
                f'Length="{step * (n - 1)!r}" Unit="m" BitInc="0" BytesInc="{inc}"/>')

    dims = dim(1, nx, dx, item) + dim(2, ny, dy, nx * item) + dim(3, nz, dz, nx * ny * nc * item)
    hw_attrs = " ".join(f"{k}={quoteattr(v)}" for k, v in img.hardware_attrs.items())
    dets = "".join("<Detector " + " ".join(f"{k}={quoteattr(v)}" for k, v in d.items()) + "/>"
                   for d in img.detectors)
    hw = (f'<Attachment Name="HardwareSetting"><ATLConfocalSettingDefinition {hw_attrs}>'
          f"<DetectorList>{dets}</DetectorList></ATLConfocalSettingDefinition></Attachment>")
    return (
        f"<Element Name={quoteattr(img.name)}><Data><Image><ImageDescription>"
        f"<Channels>{chans}</Channels><Dimensions>{dims}</Dimensions></ImageDescription>{hw}"
        f'</Image></Data><Memory Size="{img.data.nbytes}" MemoryBlockID="{block_id}"/>'
        f"<Children/></Element>"
    )


def _block(block_id: str, payload: bytes) -> bytes:
    return struct.pack("<IIBQBI", 0x70, 0, 0x2A, len(payload), 0x2A, len(block_id)) + \
        block_id.encode("utf-16-le") + payload


def write_lif(path: Path, images: list[LifTestImage], project: str = "test_project") -> Path:
    children = "".join(_image_xml(im, f"MemBlock_{i + 1}") for i, im in enumerate(images))
    xml = (
        '<LMSDataContainerHeader Version="2">'
        f'<Element Name={quoteattr(project)}><Data><Experiment/></Data>'
        '<Memory Size="0" MemoryBlockID="MemBlock_0"/>'
        f"<Children>{children}</Children></Element></LMSDataContainerHeader>"
    )
    out = struct.pack("<IIBI", 0x70, 2 * len(xml) + 5, 0x2A, len(xml)) + xml.encode("utf-16-le")
    out += _block("MemBlock_0", b"")
    for i, im in enumerate(images):
        arr = np.ascontiguousarray(im.data, dtype="<u2")
        out += _block(f"MemBlock_{i + 1}", arr.tobytes())
    path.write_bytes(out)
    return path
