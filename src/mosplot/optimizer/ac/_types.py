"""Shared data types for the reduced-MNA AC solver."""

from __future__ import annotations

from typing import NamedTuple

import numpy as np


GROUND = "gnd"
UNKNOWN_NODE = -1


class CompiledTransfer(NamedTuple):
    """Topology compiled for one voltage-transfer analysis.

    ``*_idx`` arrays hold integer node indices. ``UNKNOWN_NODE`` means the node
    is either ground or a known ideal input source, so its contribution belongs
    on the RHS instead of in the matrix. ``mos_known`` and ``cap_known`` hold
    the known-node excitation values for this analysis.
    """

    node_idx: dict[str, int]
    mos_idx: np.ndarray
    mos_known: np.ndarray
    cap_idx: np.ndarray
    cap_known: np.ndarray
    out_idx: np.ndarray
    weights: np.ndarray


class CompiledPort(NamedTuple):
    """Topology compiled for a one-port impedance query."""

    node_idx: dict[str, int]
    mos_idx: np.ndarray
    mos_known: np.ndarray
    cap_idx: np.ndarray
    cap_known: np.ndarray
    node_idx_probe: int
    reference_idx: int


class MosStamp(NamedTuple):
    """User-facing MOS stamp stored before compilation."""

    gm: float
    gmb: float
    gds: float
    d: str
    g: str
    s: str
    b: str


class CapStamp(NamedTuple):
    """User-facing capacitor stamp stored before compilation."""

    c: float
    a: str
    b: str
