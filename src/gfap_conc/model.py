"""Core data model: a validated, read-only GFAP z-stack with its provenance.

Every quantitative function in this package receives a :class:`GFAPStack`, never a
bare array, so the association between pixel values, concentration, series name,
channel and image kind (raw vs deconvolved) cannot be lost.
"""

from __future__ import annotations

import enum
import re
from dataclasses import dataclass
from typing import Any

import numpy as np


class ImageKind(str, enum.Enum):
    """Kind of image series in the .lif container."""

    RAW = "raw"
    #: ``_lng`` series: deconvolved and rescaled to the full 16-bit range (per Rose).
    #: Allowed for morphology QA only, never for intensity comparison.
    DECONVOLVED = "deconvolved"


class StackValidationError(ValueError):
    """Raised when an image array does not satisfy the (Z, Y, X) stack contract."""


class DeconvolvedIntensityError(ValueError):
    """Raised when a deconvolved (``_lng``) image is passed to an intensity analysis."""


_CONCENTRATION_RE = re.compile(r"^1:(\d+)$")


def parse_dilution(concentration: str) -> int:
    """Return the dilution factor of a label such as ``"1:500"`` (-> 500)."""
    m = _CONCENTRATION_RE.match(str(concentration).strip())
    if not m:
        raise ValueError(f"concentration must look like '1:500', got {concentration!r}")
    return int(m.group(1))


@dataclass(frozen=True, eq=False)
class GFAPStack:
    """A single-channel GFAP z-stack with explicit axes ``(Z, Y, X)``.

    The pixel array is stored as a read-only view: analysis code cannot modify
    source intensities in place. Values are kept in their native dtype; no
    rescaling of any kind is applied on construction.
    """

    series_name: str
    concentration: str
    kind: ImageKind
    channel_index: int
    data: np.ndarray
    bit_depth: int
    z_step_um: float | None = None
    pixel_size_um: tuple[float, float] | None = None  # (y, x)
    intensity_unit: str = "raw counts"
    source: str = ""

    def __post_init__(self) -> None:
        object.__setattr__(self, "kind", ImageKind(self.kind))
        parse_dilution(self.concentration)  # validates label

        data = np.asarray(self.data)
        if data.ndim != 3:
            raise StackValidationError(
                f"{self.series_name}: expected a 3-D (Z, Y, X) array, got shape {data.shape}. "
                "Select the GFAP channel and singleton T explicitly before constructing a stack."
            )
        if min(data.shape) < 1:
            raise StackValidationError(f"{self.series_name}: empty axis in shape {data.shape}")
        if not 1 <= int(self.bit_depth) <= 32:
            raise StackValidationError(f"{self.series_name}: implausible bit depth {self.bit_depth}")

        if self.kind is ImageKind.RAW:
            if not np.issubdtype(data.dtype, np.integer):
                raise StackValidationError(
                    f"{self.series_name}: raw stacks must have an integer dtype "
                    f"(got {data.dtype}); float data suggests the image was processed."
                )
            lo, hi = int(data.min()), int(data.max())
            if lo < 0:
                raise StackValidationError(f"{self.series_name}: negative raw values ({lo})")
            if hi > self.saturation_value:
                raise StackValidationError(
                    f"{self.series_name}: max value {hi} exceeds the declared {self.bit_depth}-bit "
                    f"range (max {self.saturation_value}); bit depth is wrong or data were rescaled."
                )

        if self.z_step_um is not None and not self.z_step_um > 0:
            raise StackValidationError(f"{self.series_name}: z_step_um must be > 0 (got {self.z_step_um})")

        view = data.view()
        view.flags.writeable = False
        object.__setattr__(self, "data", view)

    # ------------------------------------------------------------------ props
    @property
    def dilution_factor(self) -> int:
        return parse_dilution(self.concentration)

    @property
    def saturation_value(self) -> int:
        """Largest representable value at the *true* acquisition bit depth."""
        return 2 ** int(self.bit_depth) - 1

    @property
    def n_z(self) -> int:
        return int(self.data.shape[0])

    def z_um(self) -> np.ndarray:
        """Nominal physical depth of each slice relative to the first acquired slice."""
        if self.z_step_um is None:
            return np.full(self.n_z, np.nan)
        return np.arange(self.n_z) * float(self.z_step_um)

    def identity(self) -> dict[str, Any]:
        """Columns that tie every statistic to its source series."""
        return {
            "concentration": self.concentration,
            "dilution_factor": self.dilution_factor,
            "series_name": self.series_name,
            "image_kind": self.kind.value,
            "channel_index": self.channel_index,
            "bit_depth": int(self.bit_depth),
            "intensity_unit": self.intensity_unit,
        }


def require_raw(stack: GFAPStack, purpose: str = "intensity quantification") -> None:
    """Refuse deconvolved / full-range-normalised images for intensity work."""
    if stack.kind is not ImageKind.RAW:
        raise DeconvolvedIntensityError(
            f"{stack.series_name!r} is a {stack.kind.value} image; it was rescaled to the full "
            f"16-bit range and must not be used for {purpose}. Use the raw series instead."
        )


def sort_by_dilution(stacks: list[GFAPStack]) -> list[GFAPStack]:
    """Most concentrated (1:250) first."""
    return sorted(stacks, key=lambda s: (s.dilution_factor, s.kind.value, s.series_name))
