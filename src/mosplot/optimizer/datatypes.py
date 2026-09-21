from __future__ import annotations

import numpy as np
from dataclasses import dataclass
from typing import Callable, Literal


# A knob's role drives how the multicorner partition treats it (see evaluate.py):
#   "op"       gm/ID-like operating knob. In single-corner sizing it is chosen by the
#              optimizer; in multicorner it becomes a SOLVED unknown (re-found per corner)
#              constrained so the device width matches the frozen reference width.
#   "geom"     length / width-ratio / multiplier / passive value / mirror gain. Conserved
#              identically across every corner.
#   "external" externally applied bias held across corners (a reference current, or a bias
#              voltage from an ideal source). Optimizer-chosen in sizing, pinned in multicorner.
Role = Literal["op", "geom", "external"]


@dataclass(frozen=True)
class Spec:
    """
    Target specification for a circuit parameter.

    Attributes:
        target: The threshold value for this specification.
        mode: Optimization direction.
            "max": higher is better; penalize when actual < target (e.g. GBW, Gain, CMRR).
            "min": lower is better; penalize when actual > target (e.g. Ibias, Area).
            "eq": target an exact value; penalize deviation in both directions (e.g. VOUT_DC).
        weight: Relative importance in the cost function. Higher weight means
            violations of this spec are penalised more strongly.
        scale: Normalisation of a violation, in the spec's own unit. Needed only when the
            target is zero or negative (e.g. a saturation margin >= 0 V, where a relative
            error has no meaning); otherwise the violation is relative to |target|.
    """

    target: float
    mode: Literal["min", "max", "eq"]
    weight: float
    scale: float | None = None


@dataclass(frozen=True)
class Knob:
    """An optimizer-tuned variable in a residual circuit model.

    A circuit file declares the *physics* of a knob (name, role, which width it sets) and
    leaves `bound` as None -- the numeric search domain is PDK/spec-specific and lives in the
    config. The framework merges the two (config bounds onto circuit roles) before optimizing.

    Attributes:
        name: Variable name -- must match a name read off the State in solve_point.
        bound: Search domain (set by the config). A (lo, hi) tuple is a continuous parameter;
            a list or 1-D ndarray is a discrete choice. None in a circuit file (filled later).
        role: One of "op", "geom", "external" -- see the Role docstring. Determines how
            the multicorner partition holds or re-solves this knob.
        sets_width_of: For role="op", the device whose width this gm/ID knob sets. Used to
            build the width-match residual in multicorner. Defaults to the knob name with a
            trailing "_GMID"/"_GM_ID" stripped.
        recorner_bound: Box for this knob when a non-reference corner re-solves it (op knobs
            and RECORNER_RESOLVE knobs). The search `bound` is a design choice; at another
            corner the frozen design goes wherever the physics puts it, so it gets a wider
            box. None = the default of evaluate.recorner_bounds().
    """

    name: str
    bound: tuple[float, float] | list[float] | np.ndarray | None = None
    role: Role = "op"
    sets_width_of: str | None = None
    recorner_bound: tuple[float, float] | None = None


@dataclass(frozen=True)
class Unknown:
    """A node voltage / VDS the DC solver finds (not chosen by the optimizer).

    Attributes:
        name: Variable name -- read off the State in solve_point.
        seed: Initial guess. Either a constant, or a callable(conditions: dict) -> float
            so the seed scales with the corner (e.g. lambda c: c["vdd"] / 4) -- keep it
            PDK-agnostic by deriving from conditions, not hard-coded volts.
        bound: (lo, hi) box the solver keeps the unknown inside. Either a tuple or a
            callable(conditions: dict) -> (lo, hi) so the box can track the supply.
    """

    name: str
    seed: float | Callable[[dict], float]
    bound: tuple[float, float] | Callable[[dict], tuple[float, float]] = (0.0, 1e9)
