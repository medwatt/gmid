"""The CircuitModel contract a circuit file fills in, plus residual-scaling helpers.

A circuit subclasses `CircuitModel` and supplies six things: the topology (device /
passive / source lists + ports), the knobs, the unknowns, ``solve_point`` (all forward
lookups + matching/mirror/current substitutions, done once), ``residuals`` (the topology
equalities the solver must zero), and ``specs`` (a pure function of the solved state).
Everything else -- the solver, optimizer, corner orchestration, and netlist writer -- is
provided by the framework and never appears in a circuit file.
"""

from __future__ import annotations

from types import SimpleNamespace
from typing import ClassVar


# --- residual helpers (imported into every circuit file) ----------------------
def vres(lhs: float, rhs: float, scale: float = 0.5) -> float:
    """Voltage-equality residual, normalised so being off by `scale` volts equals 1.0."""
    return (lhs - rhs) / scale


def rres(lhs: float, rhs: float, ref: float) -> float:
    """Relative-equality residual (currents, widths), normalised by a reference magnitude."""
    return (lhs - rhs) / (abs(ref) + 1e-30)


class State(SimpleNamespace):
    """Attribute-access bag of variable values (knobs + unknowns) passed to solve_point,
    and the bag of solved quantities solve_point returns."""


class CircuitModel:
    """Base class for a residual circuit model.

    Subclasses fill the five marked sections (topology declarations, solve_point, residuals,
    specs). The framework owns solving, corners, optimizing, and netlisting.
    """

    # ---- (filled by subclass) topology / declarations ----
    NAME: ClassVar[str] = "circuit"
    PORTS: ClassVar[list] = []  # subckt port order for the netlist
    GROUND: ClassVar[str] = "vss"
    MOSFETS: ClassVar[list] = []  # list[Instance]
    PASSIVES: ClassVar[list] = []  # list[Passive]
    VSOURCES: ClassVar[list] = []  # list[VSource]
    ISOURCES: ClassVar[list] = []  # list[ISource]; netlist-only ideal current sources
    KNOBS: ClassVar[list] = []  # list[Knob]
    UNKNOWNS: ClassVar[list] = []  # list[Unknown]
    SIGNAL_NODES: ClassVar[set] = set()  # for build_ss_model (AC circuits only)
    # Knob names to ALSO re-solve in multicorner (beyond the op gm/ID knobs). Use when the
    # quantity conserved across corners is a DERIVED voltage (e.g. a tail current whose
    # applied bias is VBP = VDD - M3_VGS) rather than a directly-pinned external input.
    # Each name listed here must be pinned by a matching residual in recorner_residuals().
    RECORNER_RESOLVE: ClassVar[list] = []

    # ---- (filled by subclass) the physics ----
    def solve_point(self, v: State, dev, cond: dict) -> State:
        """Do ALL forward lookups + matching/mirror/current substitutions ONCE and return a
        State bag holding per-device DevicePoints, currents, and widths. Forward lookups
        only (dev.nmos/dev.pmos); no reverse lookup, no fixed-point loop."""
        raise NotImplementedError

    def residuals(self, b: State) -> list:
        """Return the topology equalities that must be zero at the solution (diode
        conditions, KVL/KCL closures). Use vres()/rres(); read fields off `b`."""
        raise NotImplementedError

    def specs(self, b: State, cond: dict) -> dict:
        """Pure function of the solved state: swing/margins/area/currents and, for AC
        circuits, GBW/PM/gain via build_ss_model. No solving here."""
        raise NotImplementedError

    # ---- (provided by framework; subclass overrides only as needed) ----
    def widths(self, b: State) -> dict:
        """Default: read b.W (a {device: width} dict the subclass populated)."""
        return dict(b.W)

    def lengths(self, b: State) -> dict:
        """Default: read b.L (a {device: length} dict the subclass populated)."""
        return dict(getattr(b, "L", {}))

    def device_dimensions(self, b: State) -> dict:
        """Per-device {Width, Length, Area, Current, GMID} for the report and netlist.

        Built generically from widths()/lengths() plus optional b.ID / b.GMID dicts.
        """
        W = self.widths(b)
        L = self.lengths(b)
        ID = dict(getattr(b, "ID", {}))
        GM = dict(getattr(b, "GMID", {}))
        out = {}
        for dev, w in W.items():
            length = L.get(dev)
            out[dev] = {
                "Width": w,
                "Length": length,
                "Area": (w * length) if (length is not None) else None,
                "Current": ID.get(dev),
                "GMID": GM.get(dev),
            }
        return out

    # ---- multicorner hooks (override only if RECORNER_RESOLVE is used) ----
    def freeze_extra(self, b: State) -> dict:
        """Reference values of any DERIVED quantities conserved across corners (e.g. VBP).
        Stored in frozen["extra"] and read back by recorner_residuals()."""
        return {}

    def recorner_residuals(self, b: State, frozen: dict) -> list:
        """Residuals pinning the RECORNER_RESOLVE knobs via their conserved derived quantity
        (read the reference value from frozen["extra"]). Default: none."""
        return []

    # ---- netlist hooks (override only the ones a circuit needs) ----
    def passive_values(self, b: State) -> dict:
        return {}

    def vsource_values(self, b: State, frozen: dict) -> dict:
        return {}

    def isource_values(self, b: State, frozen: dict) -> dict:
        """DC value (A) for each ISOURCES entry, keyed by name. Default: none."""
        return {}

    def mirror_currents(self, b: State) -> dict:
        """Override the auto-synthesized reference current of a mirror=-declared bias, keyed by
        the VSource name. Use when the emitted reference current is not simply ratio*Id(master)
        -- e.g. where it is the modeled IREF that makes the lossy copy deliver the
        design current. Default: none (the writer falls back to ratio*Id(master))."""
        return {}

    def netlist_context(self, corner, ref_op: State | None = None) -> dict:
        return {}

    def extra_netlist_lines(self, b: State) -> list:
        """Verbatim subcircuit lines emitted before `ends`."""
        return []
