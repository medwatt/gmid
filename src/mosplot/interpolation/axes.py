"""Axis mapping and validation for the (L, vbs, vds, gmid) grid.

All interpolation (forward and reverse) relies on np.searchsorted to bracket
query points.  searchsorted requires a strictly increasing axis -- if an axis
is unsorted, bracketing silently returns wrong indices and every lookup
produces garbage.  _require_strictly_increasing catches this at table-
construction time so the mistake is loud rather than silent.
"""

from __future__ import annotations

import numpy as np

_AXIS_DIMS = {
    "length": 0,
    "vbs": 1,
    "vds": 2,
    "gmid": 3,
}


def _require_strictly_increasing(axis, name: str) -> None:
    axis = np.asarray(axis, dtype=float)
    if axis.ndim != 1:
        raise ValueError(f"{name} axis must be 1-D.")
    if len(axis) > 1 and not np.all(np.diff(axis) > 0):
        raise ValueError(f"{name} axis must be strictly increasing.")
