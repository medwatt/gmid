from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent))

from synthetic import build_lookup_table


@pytest.fixture(scope="session")
def lookup_table() -> dict:
    return build_lookup_table()


@pytest.fixture(scope="session")
def lut_npz_path(lookup_table, tmp_path_factory) -> Path:
    """The synthetic table saved in the on-disk .npz format the loaders expect."""
    path = tmp_path_factory.mktemp("lut") / "synthetic_lut.npz"
    np.savez_compressed(path, lookup_table=np.array(lookup_table, dtype=object))
    return path
