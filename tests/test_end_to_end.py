"""End-to-end: Corner -> DeviceTable -> Optimizer -> report -> netlist file.

A 1:K NMOS current mirror sized on the synthetic square-law LUT, evaluated on
a second "slow" corner (same table with all currents scaled by 0.85), then
written out as a Spectre subcircuit. Small CMA-ES budget: the point is the
plumbing, not the quality of the optimum.
"""

import numpy as np
import pytest

from mosplot.optimizer import (
    CircuitModel,
    Corner,
    DesignReport,
    Instance,
    Knob,
    Optimizer,
    Passive,
    Spec,
    State,
    Unknown,
    VSource,
    build_ss_model,
    vres,
)
from mosplot.optimizer.netlist_writer import write_netlist

import synthetic


@pytest.fixture(scope="module")
def slow_lut_npz_path(lookup_table, tmp_path_factory):
    """A process-skewed corner: identical gm/id map, 15% less current density."""
    skewed = {}
    for name in ("nch", "pch"):
        dev = dict(lookup_table[name])
        for p in ("id", "gm", "gmbs", "gds"):
            dev[p] = dev[p] * 0.85
        skewed[name] = dev
    path = tmp_path_factory.mktemp("lut_ss") / "synthetic_lut_ss.npz"
    np.savez_compressed(path, lookup_table=np.array(skewed, dtype=object))
    return path


COND = dict(vdd=1.2, vout_dc=0.6, iref=10e-6, k=2.0, cl=1e-12)


class Mirror(CircuitModel):
    NAME = "simple_mirror"
    PORTS = ["IREF", "VOUT", "vdd", "vss"]
    GROUND = "vss"

    MOSFETS = [
        Instance("M1", "nmos", d="IREF", g="IREF", s="gnd", b="gnd"),
        Instance("M2", "nmos", d="VOUT", g="IREF", s="gnd", b="gnd"),
    ]
    PASSIVES = [Passive("CL", "cap", a="VOUT", b="gnd", external=True)]
    VSOURCES = [VSource("VDD", p="vdd", n="gnd", supply=True)]

    KNOBS = [
        Knob("M1_GMID", role="op", sets_width_of="M1"),
        Knob("M1_L", role="geom"),
    ]
    UNKNOWNS = [
        Unknown("M1_VDS", seed=lambda c: c["vdd"] / 3, bound=lambda c: (0.02, c["vdd"]))
    ]

    def solve_point(self, v, dev, cond):
        M1 = dev.nmos(gmid=v.M1_GMID, L=v.M1_L, vds=v.M1_VDS, vsb=0.0)
        M2 = dev.nmos(gmid=v.M1_GMID, L=v.M1_L, vds=cond["vout_dc"], vsb=0.0)
        ID = {"M1": cond["iref"], "M2": cond["k"] * cond["iref"]}
        W = {"M1": ID["M1"] / M1.jd, "M2": ID["M2"] / M2.jd}
        L = {"M1": v.M1_L, "M2": v.M1_L}
        GMID = {"M1": v.M1_GMID, "M2": v.M1_GMID}
        ntw = dev.nmos.table_width
        ss = {
            "M1": M1.small_signal(v.M1_GMID, ID["M1"], W["M1"], ntw, use_gmb=False),
            "M2": M2.small_signal(v.M1_GMID, ID["M2"], W["M2"], ntw, use_gmb=False),
        }
        return State(M1=M1, M2=M2, ID=ID, W=W, L=L, GMID=GMID, ss=ss, M1_VDS=v.M1_VDS)

    def residuals(self, b):
        return [vres(b.M1_VDS, b.M1.vgs)]

    def specs(self, b, cond):
        ss = build_ss_model(
            self.MOSFETS, self.PASSIVES, self.VSOURCES, b.ss, {"CL": cond["cl"]}
        )
        return {
            "Rout": ss.port("VOUT").resistance(),
            "Area": sum(b.W[d] * b.L[d] for d in b.W),
            "VDS_margin": cond["vout_dc"] - b.M2.vdsat,
        }


PARAMETERS = [
    Knob("M1_GMID", (5.0, 20.0)),
    Knob("M1_L", (150e-9, 800e-9)),
]

TARGET_SPECS = {
    "Rout": Spec(200e3, "max", 1.0),
    "Area": Spec(20e-12, "min", 0.3),
    "VDS_margin": Spec(0.1, "max", 1.0),
}


@pytest.fixture(scope="module")
def optimizer(lut_npz_path, slow_lut_npz_path):
    corners = [
        Corner("tt", str(lut_npz_path), "nch", "pch", conditions=COND, cache_dir=None),
        Corner("ss", str(slow_lut_npz_path), "nch", "pch", conditions=COND, cache_dir=None),
    ]
    opt = Optimizer(Mirror(), PARAMETERS, TARGET_SPECS, corners)
    opt.optimize(maxiter=10, n_restarts=1, seed=1)
    return opt


class TestOptimizationRun:
    def test_optimum_within_bounds(self, optimizer):
        params = optimizer.get_opt_params()
        assert 5.0 <= params["M1_GMID"] <= 20.0
        assert 150e-9 <= params["M1_L"] <= 800e-9

    def test_all_corners_solved(self, optimizer):
        assert len(optimizer.corner_results) == 2
        for r in optimizer.corner_results:
            assert r is not None
            assert r.max_residual < 1e-3
            assert set(TARGET_SPECS) <= set(r.specs)

    def test_design_is_feasible(self, optimizer):
        assert optimizer.result.fun < 1e6  # below the infeasibility floor

    def test_geometry_conserved_across_corners(self, optimizer):
        tt, ss = optimizer.corner_results
        w_tt = tt.op.W["M1"]
        w_ss = ss.op.W["M1"]
        assert w_ss == pytest.approx(w_tt, rel=1e-6)

    def test_op_knob_drifts_on_slow_corner(self, optimizer):
        # 15% less current density at the same gmid means the slow corner
        # must re-solve gm/id downward to keep the frozen width.
        gmid_ref = optimizer.get_opt_params()["M1_GMID"]
        gmid_ss = optimizer.corner_results[1].op.GMID["M1"]
        assert gmid_ss != pytest.approx(gmid_ref, rel=1e-3)

    def test_binding_corner_recorded(self, optimizer):
        assert set(optimizer.binding) <= set(TARGET_SPECS)
        assert all(b in ("tt", "ss") for b in optimizer.binding.values())


class TestReportAndNetlist:
    def test_report_contains_devices_and_specs(self, optimizer):
        text = DesignReport(optimizer).report()
        assert "M1" in text and "M2" in text
        assert "Rout" in text and "Area" in text
        assert "tt" in text and "ss" in text

    def test_netlist_written(self, optimizer, tmp_path):
        path = tmp_path / "mirror.scs"
        ref_corner = optimizer.corners[0]
        write_netlist(
            optimizer.model,
            optimizer.frozen,
            optimizer.reference_op,
            ref_corner,
            path,
        )
        text = path.read_text()
        assert "subckt simple_mirror (IREF VOUT vdd vss)" in text
        assert "M1 (IREF IREF vss vss) nch" in text
        assert "M2 (VOUT IREF vss vss) nch" in text
        assert "ends simple_mirror" in text
        assert "CL" not in text  # external load is not emitted
