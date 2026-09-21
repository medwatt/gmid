"""Write a subcircuit netlist from the solved, frozen reference design.

The geometry comes from the frozen reference operating point (the same widths/lengths the
optimizer converged on), never from a per-corner re-solve.
"""

from __future__ import annotations

from pathlib import Path

from .netlist import get_netlist_generator


def write_netlist(
    model,
    frozen,
    ref_op,
    ref_corner,
    output_path,
    *,
    netlist_format: str = "spectre",
):
    """Emit ``model``'s netlist at the frozen geometry to ``output_path``."""
    dims = {
        dev: {"Width": item["Width"], "Length": item["Length"]}
        for dev, item in frozen["dims"].items()
    }
    generator_cls = get_netlist_generator(netlist_format)
    gen = generator_cls(
        name=model.NAME,
        ports=model.PORTS,
        ground=model.GROUND,
        context=model.netlist_context(ref_corner, ref_op),
    )
    # Reference currents for auto-synthesized current-mirror biases (VSource.mirror):
    # Iref = ratio * Id(master), keyed by the bias source name so the generator can
    # emit the replica's ideal current. Merged with any explicit ISOURCES values.
    isource_params = dict(model.isource_values(ref_op, frozen))
    ref_id = dict(getattr(ref_op, "ID", {}) or {})
    mirror_curr = model.mirror_currents(ref_op)  # optional override (e.g. Idea 1's modeled IREF)
    for vs in model.VSOURCES:
        mirror = getattr(vs, "mirror", None)
        if mirror is not None:
            isource_params[vs.name] = (
                mirror_curr[vs.name]
                if vs.name in mirror_curr
                else float(getattr(vs, "ratio", 1.0)) * ref_id[mirror]
            )
    return gen.generate(
        mosfets=model.MOSFETS,
        passives=model.PASSIVES,
        vsources=model.VSOURCES,
        isources=model.ISOURCES,
        device_map={"nmos": ref_corner.nmos_name, "pmos": ref_corner.pmos_name},
        dimensions=dims,
        passive_params=model.passive_values(ref_op),
        vsource_params=model.vsource_values(ref_op, frozen),
        isource_params=isource_params,
        output_path=Path(output_path),
        extra_lines=list(model.extra_netlist_lines(ref_op)),
    )
