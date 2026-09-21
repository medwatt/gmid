from __future__ import annotations

import numpy as np

from .elements import Instance, Passive, VSource
from ..ac import SmallSignalModel, SmallSignalSolver

# Supply rails are always AC ground regardless of VSource list.
_SUPPLY_NODES: frozenset[str] = frozenset({"vdd", "vss"})

# Cache only topology, never numeric gm/gds/c values. Optimizer loops usually
# call build_ss_model thousands of times with the same node graph and changing
# device values; this cache keeps the expensive name-to-index compilation out of
# the hot path while preserving the simple builder API.
_MODEL_CACHE: dict[
    tuple[tuple[tuple[str, str, str, str], ...], tuple[tuple[str, str], ...]],
    SmallSignalModel,
] = {}


def build_ss_model(
    mosfets: list[Instance],
    passives: list[Passive],
    vsources: list[VSource],
    ss_params: dict[str, dict[str, float]],
    passive_params: dict[str, float],
    signal_nodes: set[str] | frozenset[str] | None = None,
) -> SmallSignalSolver:
    """Stamp a SmallSignalSolver from topology and computed per-element parameters.

    ac_ground is derived automatically: supply rails (vdd, vss) plus any node
    driven by a VSource, excluding nodes listed in signal_nodes. All nodes in
    ac_ground are remapped to "gnd" so the solver treats them as AC-zero.

    ss_params keys per instance: gm, gds, and optionally gmb, cgs, cgd, cdd.
    Cap values are taken as absolute values (lookup tables may return negative).
    """
    signal_nodes = frozenset(signal_nodes or ())
    ac_ground: frozenset[str] = (_SUPPLY_NODES | frozenset(vs.p for vs in vsources)) - signal_nodes
    # Netlist-only devices (e.g. a bias-mirror replica) carry no small-signal stamp.
    mosfets = [m for m in mosfets if getattr(m, "in_ss", True)]

    def _n(node: str) -> str:
        return "gnd" if node in ac_ground else node

    # Topology arrays describe where stamps connect. Value arrays carry the
    # operating-point-dependent numbers. Keeping these separate lets
    # SmallSignalModel reuse a compiled topology across candidate designs.
    mos_nodes: list[tuple[str, str, str, str]] = []
    mos_values: list[tuple[float, float, float]] = []
    cap_nodes: list[tuple[str, str]] = []
    cap_values: list[float] = []

    for inst in mosfets:
        p = ss_params[inst.name]
        mos_nodes.append((_n(inst.d), _n(inst.g), _n(inst.s), _n(inst.b)))
        mos_values.append((float(p["gm"]), float(p.get("gmb", 0.0)), float(p["gds"])))
        if "cgs" in p:
            cap_nodes.append((_n(inst.g), _n(inst.s)))
            cap_values.append(abs(float(p["cgs"])))
        if "cgd" in p:
            cap_nodes.append((_n(inst.g), _n(inst.d)))
            cap_values.append(abs(float(p["cgd"])))
        if "cdd" in p:
            cap_nodes.append((_n(inst.d), _n(inst.b)))
            cap_values.append(abs(float(p["cdd"])))

    for passive in passives:
        val = passive_params[passive.name]
        if passive.kind == "cap":
            cap_nodes.append((_n(passive.a), _n(passive.b)))
            cap_values.append(float(val))
        elif passive.kind == "res":
            mos_nodes.append((_n(passive.a), "gnd", _n(passive.b), "gnd"))
            mos_values.append((0.0, 0.0, 1.0 / float(val)))

    key = (tuple(mos_nodes), tuple(cap_nodes))
    model = _MODEL_CACHE.get(key)
    if model is None:
        model = SmallSignalModel(mos_nodes=key[0], cap_nodes=key[1])
        _MODEL_CACHE[key] = model

    return SmallSignalSolver.from_model(
        model,
        np.asarray(mos_values, dtype=float),
        np.asarray(cap_values, dtype=float),
    )


def cmfb_loop_sign(
    model,
    ref_op,
    cmfb_node: str,
    out_nodes,
    input_nodes=None,
) -> float:
    """Return the common-mode feedback loop sign (+1 / -1) for an ideal CMFB.

    The ideal CMFB drives ``cmfb_node`` to hold the average of ``out_nodes`` at a
    target. The drive must be ``vcmfb = nominal + sign * gain * (target - sensed)``
    with the right ``sign`` to be negative feedback. ``sign = sign(d V_outcm /
    d vcmfb)``, computed from the small-signal model: perturb ``cmfb_node`` (with
    the inputs held at AC ground) and read the common-mode output response.
    """
    from .elements import VSource

    pv = dict(model.passive_values(ref_op))
    for p in model.PASSIVES:
        if p.name not in pv:
            pv[p.name] = 0.0  # DC sign is independent of (open) cap values
    # Hold the differential inputs at AC ground for the CM perturbation.
    inputs = list(input_nodes if input_nodes is not None else getattr(model, "SIGNAL_NODES", ()))
    ground_inputs = [VSource(f"_cmfb_gnd{i}", p=n, n="vss") for i, n in enumerate(inputs)]
    vsources = list(model.VSOURCES) + ground_inputs
    ss = build_ss_model(
        model.MOSFETS, model.PASSIVES, vsources, ref_op.ss, pv, signal_nodes={cmfb_node}
    )
    out = {n: 1.0 / len(out_nodes) for n in out_nodes}
    g = ss.transfer(inputs={cmfb_node: 1.0}, output=out).gain()
    return 1.0 if g >= 0.0 else -1.0
