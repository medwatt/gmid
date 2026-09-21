"""Sizing, freezing, recornering and worst-casing.

One physical design is sized on a reference corner (knobs given, unknowns solved, geometry
derived) and then re-evaluated on every other corner with that geometry frozen: the op-knobs
(gm/ID) become solver unknowns again, constrained so each device's width returns to its
frozen value. Geometry and external bias stay pinned, so the system stays square and the
op-knobs are warm-started from the reference -- a handful of iterations per corner.
"""

from __future__ import annotations

from dataclasses import dataclass

from .solve import solve_point


@dataclass
class CornerResult:
    name: str
    specs: dict
    op: object = None  # the solved State
    max_residual: float = 0.0
    cost: float = float("inf")


def _op_knobs(model):
    return [k for k in model.KNOBS if k.role == "op"]


def _width_key(knob) -> str:
    return knob.sets_width_of or knob.name.replace("_GMID", "").replace("_GM_ID", "")


def _freeze(model, b, v, knob_values):
    return {
        "W": model.widths(b),
        "L": model.lengths(b),
        "dims": model.device_dimensions(b),
        "knobs": dict(knob_values),  # every knob value
        "unknown_seed": {u.name: getattr(v, u.name) for u in model.UNKNOWNS},
        "extra": model.freeze_extra(b),  # derived conserved quantities
    }


# gm/ID box for op knobs re-solved on another corner: the physical range, not the search range
RECORNER_GMID_BOUND = (1.0, 40.0)


def recorner_bounds(model) -> dict:
    """Boxes for the variables a non-reference corner re-solves.

    The config bounds describe where the optimiser may *search*; on another corner the frozen
    design can legitimately leave that box (a fixed device runs at higher gm/ID in ff, a
    margin re-derived from a fixed bias goes negative in ss). Using the search bounds there
    makes such corners look unsolvable. Defaults: op knobs get RECORNER_GMID_BOUND;
    RECORNER_RESOLVE knobs get (lo/10, 10 hi) for a positive range (currents, ratios: it never
    crosses zero), else the range widened by ten spans on each side. A knob's own
    `recorner_bound` overrides both: give one to a knob searched over a positive range that may
    change sign on another corner, such as a saturation margin.
    """
    out = {}
    for k in model.KNOBS:
        if k.role != "op" and k.name not in model.RECORNER_RESOLVE:
            continue
        rb = k.recorner_bound
        if rb is None and k.role == "op":
            rb = RECORNER_GMID_BOUND
        if rb is None and isinstance(k.bound, tuple):
            lo, hi = k.bound
            span = hi - lo
            rb = (lo / 10, hi * 10) if lo > 0 else (lo - 10 * span, hi + 10 * span)
        if rb is not None:
            out[k.name] = tuple(rb)
    return out


def check_square(model, knob_values, corners, ref_index: int = 0):
    """Raise if the sizing or the multicorner re-solve is not a square system.

    least_squares accepts an over-determined system silently and just leaves a residual, so a
    design that conserves more quantities in recorner_residuals() than it frees in
    RECORNER_RESOLVE would make every non-reference corner look infeasible. Checked once, at
    one knob point (the residual *count* does not depend on the point). Returns the counts,
    or None if the model could not be evaluated there.
    """
    try:
        res, frozen = size_reference(model, knob_values, corners[ref_index])
    except Exception:
        return None
    n_unk = len(model.UNKNOWNS)
    n_res = len(model.residuals(res.op))
    if n_res != n_unk:
        raise ValueError(f"{type(model).__name__}: residuals() returns {n_res} rows for "
                         f"{n_unk} UNKNOWNS; the reference solve is not square.")
    counts = {"sizing": (n_res, n_unk)}
    if len(corners) > 1:
        n_op = len(_op_knobs(model))
        n_cons = len(model.recorner_residuals(res.op, frozen))
        n_free = n_op + len(model.RECORNER_RESOLVE) + n_unk
        n_rows = n_res + n_op + n_cons
        if n_rows != n_free:
            raise ValueError(
                f"{type(model).__name__}: the multicorner re-solve has {n_rows} equations for "
                f"{n_free} free variables ({n_cons} conserved quantities in recorner_residuals(), "
                f"{len(model.RECORNER_RESOLVE)} knobs in RECORNER_RESOLVE). Every conserved "
                f"quantity needs one knob in RECORNER_RESOLVE that it determines."
            )
        counts["recorner"] = (n_rows, n_free)
    return counts


def size_reference(model, knob_values, ref_corner):
    """Sizing partition: knobs given, solve the UNKNOWNS, derive and freeze geometry."""
    b, v, r = solve_point(model, knob_values, ref_corner)
    frozen = _freeze(model, b, v, knob_values)
    result = CornerResult(ref_corner.name, model.specs(b, ref_corner.conditions), b, r)
    return result, frozen


def recorner(model, frozen, corner):
    """Multicorner partition: hold geometry + applied bias, re-solve the operating point.

    Op-knobs (gm/ID) are freed and constrained so each device width returns to its frozen
    value. Knobs listed in model.RECORNER_RESOLVE (e.g. a tail current whose conserved bias
    is a derived voltage) are also freed, pinned by model.recorner_residuals(). Every other
    knob (geometry, directly-applied external bias) stays at its frozen value.
    """
    op = _op_knobs(model)
    resolve = list(model.RECORNER_RESOLVE)
    free = [k.name for k in op] + resolve + [u.name for u in model.UNKNOWNS]

    knob_values = dict(frozen["knobs"])  # freed names act as placeholders, overwritten by solver

    seed = {k.name: frozen["knobs"][k.name] for k in op}
    seed.update({nm: frozen["knobs"][nm] for nm in resolve})
    seed.update(frozen["unknown_seed"])

    def extra_residuals(b, v):
        W = model.widths(b)
        r = [
            (W[_width_key(k)] - frozen["W"][_width_key(k)])
            / (abs(frozen["W"][_width_key(k)]) + 1e-30)
            for k in op
        ]
        r += list(model.recorner_residuals(b, frozen))
        return r

    b, v, r = solve_point(
        model, knob_values, corner, free_names=free, seed=seed, extra_residuals=extra_residuals,
        bounds=recorner_bounds(model),
    )
    return CornerResult(corner.name, model.specs(b, corner.conditions), b, r)


def evaluate_corners(model, knob_values, corners, ref_index: int = 0):
    """Size once on the reference corner, evaluate the frozen design on every corner.

    Returns list[CornerResult] aligned with `corners` (a slot is None if its solve failed).
    Single-corner takes the fast path: only size_reference, no width-match residuals.
    """
    try:
        ref_res, frozen = size_reference(model, knob_values, corners[ref_index])
    except Exception:
        return [None] * len(corners)

    out = [None] * len(corners)
    out[ref_index] = ref_res
    for i, c in enumerate(corners):
        if i == ref_index:
            continue
        try:
            out[i] = recorner(model, frozen, c)
        except Exception:
            out[i] = None
    return out


def worst_case(results, target_specs):
    """Per-spec worst value across corners, mode-aware, plus the binding corner.

    "max" specs take the lowest value, "min" the highest, "eq" the furthest from target.
    Non-targeted specs are passed through from the reference (first non-None) corner.

    Returns (worst, binding): worst {spec: value}, binding {spec: corner_name}.
    """
    corner_specs = [(r.name, r.specs) for r in results if r is not None and r.specs]
    worst: dict = {}
    binding: dict = {}
    for key, spec in target_specs.items():
        candidates = [(s[key], name) for name, s in corner_specs if key in s]
        if not candidates:
            continue
        if spec.mode == "max":
            val, name = min(candidates, key=lambda c: c[0])
        elif spec.mode == "min":
            val, name = max(candidates, key=lambda c: c[0])
        else:  # "eq"
            val, name = max(candidates, key=lambda c: abs(c[0] - spec.target))
        worst[key] = val
        binding[key] = name

    if corner_specs:
        for key, val in corner_specs[0][1].items():
            worst.setdefault(key, val)

    return worst, binding
