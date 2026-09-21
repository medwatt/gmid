"""Numerical kernels for reduced-MNA AC analysis."""

from __future__ import annotations

import numpy as np

from mosplot.util.numba import njit
from ._types import UNKNOWN_NODE


@njit(cache=True)
def solve_response(
    G: np.ndarray,
    C: np.ndarray,
    rhs_g: np.ndarray,
    rhs_c: np.ndarray,
    out_idx: np.ndarray,
    weights: np.ndarray,
    omega: float,
) -> complex:
    """Solve one complex frequency point and return a weighted output voltage.

    The frequency-domain system is:

        (G + j*w*C) * x = rhs_g + j*w*rhs_c

    ``rhs_c`` is non-zero when capacitors connect to known ideal input-source
    nodes. The complex matrix/RHS are filled explicitly to avoid temporary
    arrays from repeated ``astype`` calls in tight frequency-search loops.
    """

    n = G.shape[0]
    y = np.empty((n, n), dtype=np.complex128)
    rhs = np.empty(n, dtype=np.complex128)
    for i in range(n):
        rhs[i] = rhs_g[i] + 1j * omega * rhs_c[i]
        for j in range(n):
            y[i, j] = G[i, j] + 1j * omega * C[i, j]
    solution = np.linalg.solve(y, rhs)

    response = 0.0 + 0.0j
    for i in range(out_idx.shape[0]):
        response += weights[i] * solution[out_idx[i]]
    return response


@njit(cache=False)
def response_mag(
    G: np.ndarray,
    C: np.ndarray,
    rhs_g: np.ndarray,
    rhs_c: np.ndarray,
    out_idx: np.ndarray,
    weights: np.ndarray,
    omega: float,
) -> float:
    """Magnitude of the output response at one radian frequency."""

    return abs(solve_response(G, C, rhs_g, rhs_c, out_idx, weights, omega))


@njit(cache=True)
def find_unity_gain(
    G: np.ndarray,
    C: np.ndarray,
    rhs_g: np.ndarray,
    rhs_c: np.ndarray,
    out_idx: np.ndarray,
    weights: np.ndarray,
    n_iter: int,
) -> float:
    """Return the first unity-gain frequency in Hz.

    The optimizer targets compensated amplifiers, so this assumes a mostly
    low-pass response with one relevant unity crossing. The high-frequency
    bracket expansion has a large safety cap, while the bisection iteration
    count controls final crossing precision.
    """

    max_g = 1.0
    for i in range(G.shape[0]):
        for j in range(G.shape[1]):
            value = abs(G[i, j])
            if value > max_g:
                max_g = value

    min_c = np.inf
    for i in range(C.shape[0]):
        for j in range(C.shape[1]):
            value = abs(C[i, j])
            if value != 0.0 and value < min_c:
                min_c = value
    if not np.isfinite(min_c):
        min_c = 1e-12

    pole_est = max_g / min_c
    lo = 1e-3
    hi = max(pole_est * 1e3, 1.0)
    h_lo = response_mag(G, C, rhs_g, rhs_c, out_idx, weights, lo)
    if h_lo < 1.0:
        return 0.0

    h_hi = response_mag(G, C, rhs_g, rhs_c, out_idx, weights, hi)
    found_hi = False
    for _ in range(60):
        if h_hi <= 1.0:
            found_hi = True
            break
        hi *= 10.0
        h_hi = response_mag(G, C, rhs_g, rhs_c, out_idx, weights, hi)
    if not found_hi:
        return np.inf

    for _ in range(n_iter):
        mid = np.sqrt(lo * hi)
        if response_mag(G, C, rhs_g, rhs_c, out_idx, weights, mid) > 1.0:
            lo = mid
        else:
            hi = mid

    return np.sqrt(lo * hi) / (2.0 * np.pi)


@njit(cache=True)
def phase_margin(
    G: np.ndarray,
    C: np.ndarray,
    rhs_g: np.ndarray,
    rhs_c: np.ndarray,
    out_idx: np.ndarray,
    weights: np.ndarray,
    gbw_hz: float,
    dc_gain: float,
) -> float:
    """Estimate phase margin at the unity-gain frequency.

    The phase is normalized by the real DC gain to remove arbitrary
    inversion sign.
    """

    if not np.isfinite(gbw_hz) or gbw_hz <= 0.0:
        return 180.0
    if dc_gain == 0.0:
        return 180.0

    w_unity = 2.0 * np.pi * gbw_hz
    w_ref = max(w_unity * 1e-3, 1e-3)
    n_points = 1 if w_unity <= w_ref else 16
    log_ref = np.log(w_ref)
    log_unity = np.log(w_unity)

    prev_phase = 0.0
    phase_lag = 0.0
    for i in range(n_points):
        if n_points == 1:
            omega = w_unity
        else:
            omega = np.exp(log_ref + (log_unity - log_ref) * i / (n_points - 1))

        h = solve_response(G, C, rhs_g, rhs_c, out_idx, weights, omega)
        phase = np.angle(h / dc_gain)
        delta = phase - prev_phase
        while delta > np.pi:
            phase -= 2.0 * np.pi
            delta -= 2.0 * np.pi
        while delta < -np.pi:
            phase += 2.0 * np.pi
            delta += 2.0 * np.pi
        prev_phase = phase
        phase_lag = phase

    if phase_lag > 0.0:
        phase_lag -= 2.0 * np.pi

    pm = 180.0 + np.degrees(phase_lag)
    if pm < 0.0:
        return 0.0
    if pm > 180.0:
        return 180.0
    return pm


@njit(cache=True)
def stamp_admittance(
    mat: np.ndarray,
    rhs: np.ndarray,
    value: float,
    ia: int,
    ib: int,
    known_a: float,
    known_b: float,
) -> None:
    """Stamp a passive admittance between two nodes.

    Entering-current KCL is used. For row ``a`` the contribution is
    ``value * (Vb - Va)``. If either terminal is a known node, that terminal's
    voltage is moved to the RHS.
    """

    if value == 0.0:
        return

    if ia != UNKNOWN_NODE:
        mat[ia, ia] -= value
        if ib != UNKNOWN_NODE:
            mat[ia, ib] += value
        else:
            rhs[ia] -= value * known_b

    if ib != UNKNOWN_NODE:
        mat[ib, ib] -= value
        if ia != UNKNOWN_NODE:
            mat[ib, ia] += value
        else:
            rhs[ib] -= value * known_a


@njit(cache=True)
def stamp_vccs(
    G: np.ndarray,
    rhs_g: np.ndarray,
    value: float,
    op: int,
    om: int,
    cp: int,
    cm: int,
    known_cp: float,
    known_cm: float,
) -> None:
    """Stamp a current source ``value * (Vcp - Vcm)`` from ``op`` to ``om``."""

    if value == 0.0:
        return

    if op != UNKNOWN_NODE:
        if cp != UNKNOWN_NODE:
            G[op, cp] -= value
        else:
            rhs_g[op] += value * known_cp

        if cm != UNKNOWN_NODE:
            G[op, cm] += value
        else:
            rhs_g[op] -= value * known_cm

    if om != UNKNOWN_NODE:
        if cp != UNKNOWN_NODE:
            G[om, cp] += value
        else:
            rhs_g[om] -= value * known_cp

        if cm != UNKNOWN_NODE:
            G[om, cm] -= value
        else:
            rhs_g[om] += value * known_cm


@njit(cache=True)
def assemble_compiled(
    mos_values: np.ndarray,
    cap_values: np.ndarray,
    mos_idx: np.ndarray,
    mos_known: np.ndarray,
    cap_idx: np.ndarray,
    cap_known: np.ndarray,
    n_nodes: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Assemble G/C matrices and RHS vectors from a compiled topology.

    MOS value columns are ``gm, gmb, gds``. MOS node columns are ``d, g, s, b``.
    Capacitor node columns are ``a, b``. This positional convention keeps the
    hot loop compact and is documented here because it is the key extension
    point for future stamp types.
    """

    G = np.zeros((n_nodes, n_nodes), dtype=np.float64)
    C = np.zeros((n_nodes, n_nodes), dtype=np.float64)
    rhs_g = np.zeros(n_nodes, dtype=np.float64)
    rhs_c = np.zeros(n_nodes, dtype=np.float64)

    for i in range(mos_values.shape[0]):
        gm = mos_values[i, 0]
        gmb = mos_values[i, 1]
        gds = mos_values[i, 2]
        d = mos_idx[i, 0]
        g = mos_idx[i, 1]
        s = mos_idx[i, 2]
        b = mos_idx[i, 3]

        stamp_admittance(
            G,
            rhs_g,
            gds,
            d,
            s,
            mos_known[i, 0],
            mos_known[i, 2],
        )
        stamp_vccs(
            G,
            rhs_g,
            gm,
            d,
            s,
            g,
            s,
            mos_known[i, 1],
            mos_known[i, 2],
        )
        stamp_vccs(
            G,
            rhs_g,
            gmb,
            d,
            s,
            b,
            s,
            mos_known[i, 3],
            mos_known[i, 2],
        )

    for i in range(cap_values.shape[0]):
        stamp_admittance(
            C,
            rhs_c,
            cap_values[i],
            cap_idx[i, 0],
            cap_idx[i, 1],
            cap_known[i, 0],
            cap_known[i, 1],
        )

    return G, C, rhs_g, rhs_c
