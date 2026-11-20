"""Synthetic GFAP-like stacks for tests and the pipeline demo.

These are NOT models of real tissue. They only provide known ground truth
(signal scale, background level, a mid-depth intensity dip) so that the code
path can be checked end to end without lab data. All figures generated from them
are watermarked "SYNTHETIC".
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi
from skimage.draw import line

from .model import GFAPStack, ImageKind


def process_skeleton(
    shape: tuple[int, int, int], n_processes: int, rng: np.random.Generator, full_depth: bool = True
) -> np.ndarray:
    """Random straight 'processes', before blurring.

    With ``full_depth`` every process spans all slices, so any depth trend comes only
    from :func:`depth_attenuation`. Without it, processes occupy random z-ranges -
    a reminder that uneven structure density alone can create a depth profile.
    """
    nz, ny, nx = shape
    vol = np.zeros(shape, dtype=np.float64)
    for _ in range(n_processes):
        r0, c0 = rng.integers(0, ny), rng.integers(0, nx)
        r1, c1 = rng.integers(0, ny), rng.integers(0, nx)
        rr, cc = line(int(r0), int(c0), int(r1), int(c1))
        if full_depth:
            za, zb = 0, nz
        else:
            za = int(rng.integers(0, nz))
            zb = int(min(nz, za + rng.integers(max(2, nz // 4), nz + 1)))
        weight = rng.uniform(0.4, 1.0)
        vol[za:zb, rr, cc] = np.maximum(vol[za:zb, rr, cc], weight)
    return vol


def depth_attenuation(nz: int, middle_dip: float = 0.4, first_side_boost: float = 0.0) -> np.ndarray:
    """Multiplicative per-slice factor: Gaussian dip of depth ``middle_dip`` at mid-stack,
    optionally brighter near the first acquired slice."""
    z = np.arange(nz)
    mid = (nz - 1) / 2
    width = max(nz / 5, 1.0)
    f = 1.0 - middle_dip * np.exp(-0.5 * ((z - mid) / width) ** 2)
    if first_side_boost:
        f *= 1.0 + first_side_boost * np.exp(-z / max(nz / 6, 1.0))
    return f


def make_stack_array(
    shape: tuple[int, int, int] = (24, 128, 128),
    signal_amplitude: float = 1500.0,
    background_level: float = 60.0,
    n_processes: int = 40,
    middle_dip: float = 0.4,
    first_side_boost: float = 0.0,
    bit_depth: int = 12,
    structure_seed: int = 1,
    noise_seed: int = 2,
    poisson: bool = True,
) -> np.ndarray:
    """Integer stack = background + amplitude * blurred processes * depth factor (+ Poisson noise)."""
    struct = process_skeleton(shape, n_processes, np.random.default_rng(structure_seed))
    struct = ndi.gaussian_filter(struct, sigma=(0.7, 1.2, 1.2))
    if struct.max() > 0:
        struct /= struct.max()
    expected = background_level + signal_amplitude * struct * depth_attenuation(shape[0], middle_dip, first_side_boost)[:, None, None]
    rng = np.random.default_rng(noise_seed)
    values = rng.poisson(expected) if poisson else np.rint(expected)
    return np.clip(values, 0, 2**bit_depth - 1).astype(np.uint16)


def make_concentration_series(
    amplitudes: dict[str, float] | None = None,
    shape: tuple[int, int, int] = (24, 128, 128),
    bit_depth: int = 12,
    background_level: float = 60.0,
    include_deconvolved: bool = True,
    z_step_um: float = 0.5,
    pixel_size_um: float = 0.18,
) -> list[GFAPStack]:
    """Four raw synthetic stacks (shared structure, different amplitude) plus optional
    ``_lng``-style deconvolved copies rescaled to the full 16-bit range."""
    amplitudes = amplitudes or {"1:250": 2000.0, "1:500": 1600.0, "1:1000": 1100.0, "1:2000": 600.0}
    out = []
    for i, (conc, amp) in enumerate(amplitudes.items()):
        arr = make_stack_array(shape, amp, background_level, bit_depth=bit_depth, noise_seed=10 + i)
        name = f"SYNTH_GFAP_{conc.replace(':', '-')}"
        out.append(
            GFAPStack(name, conc, ImageKind.RAW, 0, arr, bit_depth, z_step_um, (pixel_size_um, pixel_size_um),
                      "raw counts (synthetic)", "synthetic")
        )
        if include_deconvolved:
            f = arr.astype(np.float64)
            rescaled = np.rint((f - f.min()) / max(f.max() - f.min(), 1) * 65535).astype(np.uint16)
            out.append(
                GFAPStack(f"{name}_lng", conc, ImageKind.DECONVOLVED, 0, rescaled, 16, z_step_um,
                          (pixel_size_um, pixel_size_um), "rescaled (synthetic)", "synthetic")
            )
    return out
