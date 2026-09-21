"""The run() entry point a config file calls.

Instantiate the circuit model, optimize it across the corners, print the report, and write
the netlist for the frozen reference design.
"""

from __future__ import annotations

from .design_report import DesignReport
from .netlist_writer import write_netlist
from .optimizer import Optimizer


def run(
    *,
    circuit,
    parameters=None,
    target_specs,
    corners,
    output_module,
    maxiter: int = 300,
    seed: int = 1,
    n_restarts: int = 1,
    ref_index: int = 0,
    netlist_format: str = "spectre",
    workers: int = 1,
):
    model = circuit()
    opt = Optimizer(
        model,
        parameters or model.KNOBS,
        target_specs,
        corners,
        ref_index=ref_index,
    )
    opt.optimize(maxiter=maxiter, seed=seed, n_restarts=n_restarts, workers=workers)
    print(DesignReport(opt).report())
    write_netlist(
        model,
        opt.frozen,
        opt.reference_op,
        corners[ref_index],
        output_module,
        netlist_format=netlist_format,
    )
    return opt
