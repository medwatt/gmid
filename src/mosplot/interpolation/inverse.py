"""Pure 1-D inverse solver: find x where y(x) == target.

No MOSFET logic. No Expression logic. Just 1D array inversion.

Key design decisions:
  - Full curve scan, not bisection -- the default solver must handle
    non-monotone curves (e.g. gm, gds, fT, gm/gds vs gmid).
  - NaN gaps split the curve into independent finite segments.
    Crossings are never interpolated across non-finite gaps.
  - Multiple crossings are detected, warned about, and resolved
    by mode ("first", "last", "nearest", "all").

The optimizer hot path uses _inverse_1d_nb (Numba JIT in interp.py)
for the common case (monotone, first crossing, clip); this module
handles the remaining modes, warnings, and edge cases.
"""

from __future__ import annotations

import warnings
import numpy as np

_VALID_MODES = {"first", "last", "nearest", "all"}
_VALID_OUT_OF_RANGE = {"nan", "clip", "raise"}


def _finite_segments(valid: np.ndarray) -> list[tuple[int, int]]:
    """Return half-open [start, stop) ranges of contiguous True values.

    A single NaN/inf breaks the curve.  A segment can be length 1
    (exact grid hit) though sign-change interpolation needs ≥ 2 points.
    """
    segments = []
    n = len(valid)
    i = 0

    while i < n:
        while i < n and not valid[i]:
            i += 1
        start = i
        while i < n and valid[i]:
            i += 1
        stop = i
        if stop > start:
            segments.append((start, stop))

    return segments


def _dedupe_sorted(values: list[float], *, atol: float, rtol: float) -> list[float]:
    """Sort and remove near-duplicate crossing positions.

    Needed because a plateau of exact hits and an adjacent sign-change
    crossing can produce the same interpolated x value up to round-off.
    """
    if not values:
        return []
    values = sorted(values)
    out = [values[0]]
    for v in values[1:]:
        if not np.isclose(v, out[-1], atol=atol, rtol=rtol):
            out.append(v)
    return out


def _inverse_1d(
    x,
    y,
    target,
    *,
    mode: str = "first",
    out_of_range: str = "nan",
    warn_multiple: bool = True,
    rtol: float = 1e-9,
    atol: float = 0.0,
) -> float | list[float]:
    """Find x where y(x) == target by linear crossing scan.

    For each contiguous finite segment:
      1. Collect exact grid-point hits (np.isclose).
      2. Find sign-change crossings between adjacent non-exact points.
      3. Linearly interpolate each crossing.

    After scanning all segments:
      - Deduplicate near-identical solutions.
      - If multiple crossings remain and mode != "all", issue a warning.
      - Return according to mode.

    Parameters
    ----------
    x, y : array-like
        Same-length 1-D arrays.  Non-finite values break the curve.
    target : float
    mode : {"first", "last", "nearest", "all"}
    out_of_range : {"nan", "clip", "raise"}
    warn_multiple : bool
    rtol, atol : float
        Tolerances for exact-hit detection.

    Returns
    -------
    float or list[float]
    """
    if mode not in _VALID_MODES:
        raise ValueError(
            f"Unknown inverse lookup mode: {mode!r}. "
            f"Expected one of {sorted(_VALID_MODES)}."
        )

    if out_of_range not in _VALID_OUT_OF_RANGE:
        raise ValueError(
            f"Unknown out_of_range policy: {out_of_range!r}. "
            f"Expected one of {sorted(_VALID_OUT_OF_RANGE)}."
        )

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    target = float(target)

    if x.shape != y.shape:
        raise ValueError(f"x and y must have the same shape, got {x.shape} and {y.shape}.")

    if x.ndim != 1:
        raise ValueError(f"x and y must be 1-D arrays, got ndim={x.ndim}.")

    valid = np.isfinite(x) & np.isfinite(y)

    if not np.any(valid):
        if out_of_range == "raise":
            raise ValueError("Reverse lookup has no finite samples.")
        return [] if mode == "all" else np.nan

    solutions = []

    # ── scan each contiguous finite segment independently ──
    for start, stop in _finite_segments(valid):
        xs = x[start:stop]
        ys = y[start:stop]

        if xs.size == 0:
            continue

        diff = ys - target
        close = np.isclose(ys, target, rtol=rtol, atol=atol)

        # Exact grid-point hits (including plateau points).
        for idx in np.where(close)[0]:
            solutions.append(float(xs[idx]))

        # Sign-change crossings between adjacent non-exact finite points.
        for i in range(xs.size - 1):
            if close[i] or close[i + 1]:
                continue

            d0 = diff[i]
            d1 = diff[i + 1]

            if d0 * d1 < 0.0:                     # opposite signs -> crossing
                y0 = ys[i]
                y1 = ys[i + 1]

                if y1 != y0:                      # linear interpolation
                    t = (target - y0) / (y1 - y0)
                    solutions.append(float(xs[i] + t * (xs[i + 1] - xs[i])))

    solutions = _dedupe_sorted(solutions, atol=atol, rtol=rtol)

    # ── handle no-solution cases ──
    finite_y = y[valid]
    finite_x = x[valid]
    finite_diff = finite_y - target

    if not solutions:
        if out_of_range == "clip":
            return float(finite_x[int(np.nanargmin(np.abs(finite_diff)))])
        if out_of_range == "raise":
            raise ValueError(
                f"Target {target} is outside the achievable finite range "
                f"[{float(np.nanmin(finite_y)):.4g}, {float(np.nanmax(finite_y)):.4g}]."
            )
        return [] if mode == "all" else np.nan

    # ── multiple-crossing warning ──
    if len(solutions) > 1 and mode != "all" and warn_multiple:
        warnings.warn(
            f"Reverse lookup found {len(solutions)} crossings: "
            f"{[f'{s:.4g}' for s in solutions]}. "
            f"Returning the {mode!r} crossing. "
            f"Pass mode='all' to retrieve all solutions, or pass "
            f"warn_multiple=False if this ambiguity is expected.",
            UserWarning,
            stacklevel=3,
        )

    # ── resolve by mode ──
    if mode == "all":
        return solutions

    if mode == "first":
        return solutions[0]

    if mode == "last":
        return solutions[-1]

    if mode == "nearest":
        nearest_idx = int(np.nanargmin(np.abs(finite_diff)))
        x_nearest = float(finite_x[nearest_idx])
        return min(solutions, key=lambda s: abs(s - x_nearest))

    raise RuntimeError("unreachable")
