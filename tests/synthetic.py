"""Synthetic square-law MOSFET lookup table for tests.

Builds a lookup table in the exact layout the simulator backends produce
(data shaped (n_L, n_vbs, n_vgs, n_vds)), but from a closed-form square-law
model so every interpolated value can be checked against an analytic truth:

    vth(vbs) = VTH0 + GAMMA * (sqrt(PHI - vbs) - sqrt(PHI))
    vov      = |vgs| - |vth|
    id       = 0.5 * KP * (W / L) * vov^2 * (1 + LAMBDA(L) * |vds|)
    gm       = KP * (W / L) * vov * (1 + LAMBDA(L) * |vds|)
    gm/id    = 2 / vov

The saturation equation is used at every bias point (no triode branch); the
tests only need a smooth, invertible surface, not physical accuracy.

PMOS entries follow the on-disk conventions the optimizer's DeviceTable
expects: vgs/vds/vth/vdsat stored negative, id/gm stored positive. The vbs axis is
at or below zero for both polarities here; build_device_vbs_up() gives the PMOS
convention of the LUT generator (vbs at or above zero).
"""

from __future__ import annotations

import numpy as np

WIDTH = 10e-6
KP = 200e-6
VTH0 = 0.4
GAMMA = 0.3
PHI = 0.8
LAMBDA0 = 0.06  # at L = 400 nm; scales as 1/L

LENGTHS = np.array([100e-9, 200e-9, 400e-9, 800e-9])
VBS_AXIS = np.array([-0.6, -0.3, 0.0])
VGS_AXIS = np.linspace(0.0, 1.2, 61)
VDS_AXIS = np.linspace(0.0, 1.2, 21)

PARAMETER_NAMES = ["id", "vth", "vdsat", "gm", "gmbs", "gds", "cgs", "cgd", "cdd"]


def lam(length: float | np.ndarray) -> float | np.ndarray:
    return LAMBDA0 * 400e-9 / length


def vth_of_vbs(vbs: float | np.ndarray) -> float | np.ndarray:
    return VTH0 + GAMMA * (np.sqrt(PHI - vbs) - np.sqrt(PHI))


def analytic_point(length: float, vgs: float, vds: float, vbs: float) -> dict:
    """Magnitude-domain analytic truth at one bias point."""
    vth = vth_of_vbs(vbs)
    vov = max(vgs - vth, 0.0)
    klw = KP * WIDTH / length
    boost = 1.0 + lam(length) * vds
    id_ = 0.5 * klw * vov**2 * boost
    gm = klw * vov * boost
    cox_area = 5e-3 * WIDTH * length
    return {
        "vth": vth,
        "vov": vov,
        "id": id_,
        "gm": gm,
        "gmbs": 0.2 * gm,
        "gds": 0.5 * klw * vov**2 * lam(length),
        "vdsat": vov,
        "gmid": (2.0 / vov) if vov > 0 else np.nan,
        "cgs": (2.0 / 3.0) * cox_area,
        "cgd": 0.1 * cox_area,
        "cdd": 0.2 * cox_area,
    }


def analytic_from_gmid(length: float, gmid: float, vds: float, vbs: float) -> dict:
    """Analytic truth parameterised by gm/Id instead of vgs."""
    vov = 2.0 / gmid
    vgs = vov + vth_of_vbs(vbs)
    return analytic_point(length, vgs, vds, vbs)


def _device_arrays(polarity: str) -> dict:
    """Dense (n_L, n_vbs, n_vgs, n_vds) parameter arrays."""
    L = LENGTHS[:, None, None, None]
    vbs = VBS_AXIS[None, :, None, None]
    vgs = VGS_AXIS[None, None, :, None]
    vds = VDS_AXIS[None, None, None, :]

    vth = vth_of_vbs(vbs)
    vov = np.maximum(vgs - vth, 0.0)
    klw = KP * WIDTH / L
    boost = 1.0 + lam(L) * vds
    id_ = 0.5 * klw * vov**2 * boost
    gm = klw * vov * boost
    gds = 0.5 * klw * vov**2 * lam(L)
    cox_area = 5e-3 * WIDTH * L
    shape = np.broadcast_shapes(id_.shape, (1, 1, 1, 1))

    sign = -1.0 if polarity == "p" else 1.0
    return {
        # id/gm stay positive for both polarities so gm/id > 0.
        "id": np.broadcast_to(id_, shape).copy(),
        "gm": np.broadcast_to(gm, shape).copy(),
        "gmbs": np.broadcast_to(0.2 * gm, shape).copy(),
        "gds": np.broadcast_to(gds, shape).copy(),
        # Voltage-like quantities carry the polarity sign on disk.
        "vth": np.broadcast_to(sign * vth, shape).copy(),
        "vdsat": np.broadcast_to(sign * vov, shape).copy(),
        "cgs": np.broadcast_to(cox_area, shape).copy(),
        "cgd": np.broadcast_to(0.1 * cox_area, shape).copy(),
        "cdd": np.broadcast_to(0.2 * cox_area, shape).copy(),
    }


def build_device(polarity: str) -> dict:
    sign = -1.0 if polarity == "p" else 1.0
    dev = _device_arrays(polarity)
    dev["length"] = LENGTHS.copy()
    dev["vbs"] = VBS_AXIS.copy()
    dev["vgs"] = sign * VGS_AXIS
    dev["vds"] = sign * VDS_AXIS
    dev["parameter_names"] = list(PARAMETER_NAMES)
    dev["device_parameters"] = {"w": WIDTH}
    dev["model_name"] = "nch" if polarity == "n" else "pch"
    return dev


def build_device_vbs_up(polarity: str = "p") -> dict:
    """Same device with the vbs axis the LookupTableGenerator writes for PMOS: at or above
    zero (reverse bias raises the bulk above the source), increasing. Data reordered to match."""
    dev = build_device(polarity)
    for k, v in list(dev.items()):
        if isinstance(v, np.ndarray) and v.ndim == 4:
            dev[k] = v[:, ::-1].copy()
    dev["vbs"] = (-VBS_AXIS)[::-1].copy()
    return dev


def build_lookup_table() -> dict:
    return {
        "nch": build_device("n"),
        "pch": build_device("p"),
        "description": "synthetic square-law test table",
        "simulator": "synthetic",
        "parameter_names": list(PARAMETER_NAMES),
        "device_parameters": {"w": WIDTH},
    }
