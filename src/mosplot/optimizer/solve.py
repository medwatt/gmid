"""The generic DC solver: drive model.residuals(...) == 0 for the free unknowns.

The solver is deliberately topology-agnostic -- it never knows what the equations mean. The
same routine runs single-corner sizing (free = the unknowns, knobs pinned) and multicorner
recornering (free = op-knobs + unknowns, geometry pinned), differing only in which variables
are handed to it as free and which extra residuals are appended.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import least_squares

from .template import State


def solve_point(
    model,
    knob_values: dict,
    corner,
    *,
    free_names=None,
    pinned=None,
    seed=None,
    extra_residuals=None,
    bounds=None,
    max_nfev: int = 2000,
):
    """Solve model.residuals(...) == 0 for the free unknowns on a given corner.

    Args:
        model: the CircuitModel instance.
        knob_values: {name: value} for all knobs (some may be in free_names or pinned).
        corner: a Corner; exposes .nmos/.pmos DeviceTables and .conditions.
        free_names: names solved by least_squares. Default: every Unknown name.
        pinned: {name: value} forced fixed (overrides knob_values / unknown seeds).
        seed: {name: value} warm-start for free_names (else use the Unknown seed / knob value).
        extra_residuals: optional callable(b, v) -> list[float] appended to model.residuals(b)
            (used for the multicorner width-match constraints).
        bounds: optional {name: (lo, hi)} overriding the box of those free variables (the
            multicorner re-solve uses it to replace the knobs' search bounds).
        max_nfev: least_squares iteration cap.

    Returns:
        (b, v, max_abs_residual): b = solve_point() output, v = State of all values.
    """
    cond = dict(corner.conditions)
    pinned = pinned or {}
    seed = seed or {}
    if free_names is None:
        free_names = [u.name for u in model.UNKNOWNS]

    if not free_names:
        v = State(**knob_values)
        for nm, val in pinned.items():
            setattr(v, nm, val)
        b = model.solve_point(v, corner, cond)
        return b, v, 0.0

    bounds_by_name = {
        u.name: tuple(u.bound(cond) if callable(u.bound) else u.bound) for u in model.UNKNOWNS
    }
    bounds_by_name.update(
        {k.name: tuple(k.bound) for k in model.KNOBS if isinstance(k.bound, tuple)}
    )
    bounds_by_name.update(bounds or {})

    x0, lo, hi = [], [], []
    for nm in free_names:
        if nm in seed:
            guess = seed[nm]
        else:
            u = next((u for u in model.UNKNOWNS if u.name == nm), None)
            if u is not None:
                guess = u.seed(cond) if callable(u.seed) else u.seed
            else:
                guess = knob_values.get(nm, 1.0)
        b_lo, b_hi = bounds_by_name.get(nm, (-1e9, 1e9))
        x0.append(min(max(guess, b_lo), b_hi))
        lo.append(b_lo)
        hi.append(b_hi)

    def assemble(xfree) -> State:
        v = State(**knob_values)
        for nm, val in zip(free_names, xfree):
            setattr(v, nm, val)
        for nm, val in pinned.items():
            setattr(v, nm, val)
        return v

    def fn(xfree):
        v = assemble(xfree)
        b = model.solve_point(v, corner, cond)
        r = list(model.residuals(b))
        if extra_residuals is not None:
            r += list(extra_residuals(b, v))
        return np.nan_to_num(np.asarray(r, dtype=float), nan=1e3, posinf=1e3, neginf=1e3)

    res = least_squares(
        fn,
        x0,
        bounds=(lo, hi),
        method="trf",
        x_scale="jac",
        ftol=1e-12,
        xtol=1e-12,
        max_nfev=max_nfev,
    )
    v = assemble(res.x)
    b = model.solve_point(v, corner, cond)
    return b, v, float(np.max(np.abs(res.fun)))
