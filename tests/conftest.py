import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from kubeviz.synth import make_synthetic_cube  # noqa: E402


@pytest.fixture(scope="session")
def synth_cube(tmp_path_factory):
    d = tmp_path_factory.mktemp("data")
    fname = str(d / "synth_cube.fits")
    _, truth = make_synthetic_cube(fname, nx=12, ny=10, redshift=0.02, nan_corner=True)
    return fname, truth
