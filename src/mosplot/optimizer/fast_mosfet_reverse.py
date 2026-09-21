"""Expression-based reverse lookup helpers for FastMosfet.

These connect GmIdTable.lookup_curve (raw-parameter interpolation along one
free axis) with Expression evaluation and 1-D inversion to answer:

  "Find one axis value such that expression(...) == target."

Flow for each reverse lookup:
  1. lookup_curve:  extract the full 1-D curve for each raw parameter
                    needed by the expression (e.g. "id" for JD).
  2. _expression_curve:  evaluate the Expression over the 1-D curves,
                         producing a single y(axis) array.
  3. inverse solver: find crossing(s) of y(axis) == target.

The optimizer's hot path (mode="first", out_of_range="clip", no warning)
takes the Numba JIT fast path (_inverse_1d_nb) for ~26 µs/call.
All other modes/options fall back to the Python _inverse_1d in inverse.py.
"""

from __future__ import annotations

import numpy as np

from mosplot.interpolation.inverse import _inverse_1d
from mosplot.interpolation.fast.interp import _inverse_1d_nb

# The four grid axes.  Used to distinguish raw table parameters from
# axis-level variables in Expression evaluation.
_AXIS_NAMES = frozenset({"length", "gmid", "vds", "vbs"})


def _expression_curve(
    *,
    solve_for: str,
    expression,
    axis: np.ndarray,
    raw_curves: dict[str, np.ndarray],
    length: float | None,
    gmid: float | None,
    vds: float | None,
    vbs: float | None,
) -> np.ndarray:
    """Evaluate one Expression over a 1-D curve along the solved axis.

    Each variable in the expression is resolved in order:
      1. If it's a raw table parameter (e.g. "id", "gm") -> raw_curves[var]
         (a 1-D NumPy array over the axis points).
      2. If it's the solved axis itself -> axis (the x-values).
      3. If it's a fixed axis -> a constant array np.full_like(axis, value).

    The expression must support vectorized NumPy inputs -- an Expression
    designed for scalar forward lookup will raise ValueError here.
    """
    fixed_values = {
        "length": length,
        "gmid": gmid,
        "vds": vds,
        "vbs": vbs,
    }

    args = []

    for var in expression.variables:
        if var in raw_curves:  # raw table param
            args.append(raw_curves[var])
        elif var == solve_for:  # the free axis
            args.append(axis)
        elif var in fixed_values and fixed_values[var] is not None:  # fixed axis
            args.append(np.full_like(axis, fixed_values[var], dtype=float))
        else:
            raise KeyError(
                f"Expression variable {var!r} is not available as a raw "
                f"parameter, solved axis, or fixed axis."
            )

    # Let NumPy errors (TypeError, etc.) propagate with a clear message.
    try:
        y = expression.function(*args)
    except Exception:
        raise ValueError(
            f"Expression {expression!r} could not be evaluated over a NumPy "
            f"array curve. Expressions used in reverse lookup must support "
            f"vectorized NumPy array inputs."
        ) from None

    y = np.asarray(y, dtype=float)

    # Scalar result -> broadcast to match the axis.
    if y.shape == ():
        return np.full_like(axis, float(y), dtype=float)

    if y.shape != axis.shape:
        raise ValueError(
            f"Expression returned shape {y.shape}, expected {axis.shape}. "
            f"Expressions used in reverse lookup must support NumPy array input."
        )

    return y


def _lookup_axis_for_impl(
    self,
    solve_for: str,
    *,
    expression,
    target: float,
    length: float | None = None,
    gmid: float | None = None,
    vds: float | None = None,
    vbs: float | None = None,
    mode: str = "first",
    out_of_range: str = "nan",
    warn_multiple: bool = True,
    method: str = "scan",
) -> float | list[float]:
    """Core reverse-lookup implementation (separated from FastMosfet for testing).

    Steps:
      1. Validate inputs (single expression, known axis, solved axis omitted).
      2. Collect raw table parameters needed by the expression.
      3. Extract 1-D curves via GmIdTable.lookup_curve.
      4. Evaluate the expression over the curve -> y(axis).
      5. Solve y(axis) == target via inverse solver.

    Two solvers are available (see FastMosfet.lookup_axis_for docstring):
      method="scan" -> Python _inverse_1d (full-featured, all modes/warnings)
      method="fast" -> Numba _inverse_1d_nb (~4× faster, restricted options)
    """

    # ── validate inputs ──
    if isinstance(expression, (list, tuple)):
        raise TypeError(
            "lookup_axis_for() expects a single Expression, not a list. "
            "Reverse lookup solves one equation at a time."
        )

    if solve_for not in _AXIS_NAMES:
        raise ValueError(
            f"Unknown solve axis: {solve_for!r}. Must be one of {sorted(_AXIS_NAMES)}."
        )

    if method not in ("scan", "fast"):
        raise ValueError(f"Unknown method: {method!r}.  Use 'scan' (default) or 'fast'.")

    values = {
        "length": length,
        "gmid": gmid,
        "vds": vds,
        "vbs": vbs,
    }

    if values[solve_for] is not None:
        raise ValueError(f"Axis {solve_for!r} is being solved for and must not be provided.")

    # ── validate fast-mode restrictions ──
    if method == "fast":
        if mode != "first":
            raise ValueError(
                f"method='fast' requires mode='first', got mode={mode!r}. "
                f"Use method='scan' for mode='last'/'nearest'/'all'."
            )
        if out_of_range not in ("clip", "nan"):
            raise ValueError(
                f"method='fast' requires out_of_range in ('clip', 'nan'), "
                f"got {out_of_range!r}.  Use method='scan' for 'raise'."
            )

    # ── extract 1-D curves for raw parameters needed by the expression ──
    needed = sorted(set(expression.variables) - _AXIS_NAMES)

    axis, raw_curves = self._table.lookup_curve(
        solve_for=solve_for,
        params=needed,
        length=length,
        gmid=gmid,
        vds=vds,
        vbs=vbs,
    )

    # ── evaluate the expression along the curve ──
    y = _expression_curve(
        solve_for=solve_for,
        expression=expression,
        axis=axis,
        raw_curves=raw_curves,
        length=length,
        gmid=gmid,
        vds=vds,
        vbs=vbs,
    )

    # ── solve y(axis) == target ──
    if method == "fast":
        oor = 1 if out_of_range == "clip" else 0
        return _inverse_1d_nb(axis, y, target, oor)

    return _inverse_1d(
        axis,
        y,
        target,
        mode=mode,
        out_of_range=out_of_range,
        warn_multiple=warn_multiple,
    )
