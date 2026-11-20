import numpy as np
import pytest

from gfap_conc.model import GFAPStack, ImageKind


def make_stack(data, concentration="1:500", kind=ImageKind.RAW, bit_depth=12, z_step_um=1.0,
               name=None, pixel_size_um=(0.2, 0.2)) -> GFAPStack:
    return GFAPStack(
        series_name=name or f"S_{concentration}_{kind.value if isinstance(kind, ImageKind) else kind}",
        concentration=concentration,
        kind=kind,
        channel_index=0,
        data=np.asarray(data),
        bit_depth=bit_depth,
        z_step_um=z_step_um,
        pixel_size_um=pixel_size_um,
    )


@pytest.fixture
def rng():
    return np.random.default_rng(1234)
