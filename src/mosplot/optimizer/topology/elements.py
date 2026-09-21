from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class Instance:
    """One transistor instance in the circuit topology.

    Set in_ss=False for devices that belong to the generated netlist but NOT to
    the small-signal model -- e.g. a bias-mirror replica whose only job is to set
    a DC gate voltage. Such a device is sized/area-counted and emitted, but the
    AC solver ignores it (the bias node is still AC-grounded by its VSource
    stand-in).
    """

    name: str
    kind: str  # device type key ("nmos" or "pmos")
    d: str
    g: str
    s: str
    b: str
    in_ss: bool = field(default=True)


@dataclass(frozen=True)
class Passive:
    """A two-terminal passive element (capacitor or resistor).

    Set external=True for elements that belong to the AC model but not to the
    generated subcircuit netlist (e.g. an off-chip load capacitor).
    """

    name: str
    kind: str  # "cap" | "res"
    a: str
    b: str
    external: bool = field(default=False)


@dataclass(frozen=True)
class VSource:
    """An ideal DC voltage source inside the subcircuit.

    The positive terminal p is treated as AC ground by the small-signal solver
    because it is driven by an ideal source.

    Set emit=False to keep the node AC-grounded in the small-signal model while
    NOT writing the source into the netlist -- used when a bias node is realized
    by some other network (e.g. a current-mirror replica) in the generated
    circuit but is just an ideal voltage stand-in for the AC analysis.

    Set mirror="<device>" to declare that this bias node is the gate of a
    current-source device and is realized by a current mirror. The netlist
    generator then auto-synthesizes a diode-connected replica of that device fed
    by an ideal reference current (the device's own Id, scaled by `ratio`), so
    the bias node self-biases across corners exactly like a real mirror -- and
    the source is NOT written as an ideal voltage. The small-signal model still
    treats p as AC ground (the replica is a DC bias network).
    """

    name: str
    p: str  # driven (bias) node
    n: str  # reference node (typically "vss")
    supply: bool = field(default=False)
    emit: bool = field(default=True)
    mirror: str | None = field(default=None)  # device name this bias mirrors
    ratio: float = field(default=1.0)  # reference-current scale (Iref = ratio * Id_master)


@dataclass(frozen=True)
class ISource:
    """An ideal DC current source inside the subcircuit.

    Netlist-only: it has no small-signal stamp (an ideal current source is an AC
    open). Positive dc sinks current from p to n, so to pull a reference current
    out of a bias node use p=<bias node>, n="vss".
    """

    name: str
    p: str
    n: str
    supply: bool = field(default=False)
