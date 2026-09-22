from __future__ import annotations

from mosplot.optimizer import (
    CircuitModel,
    Instance,
    Knob,
    Passive,
    State,
    Unknown,
    VSource,
    run,
    vres,
)


class Circuit(CircuitModel):
    NAME = "simple_mirror"
    PORTS = ["IREF", "VOUT", "vdd", "vss"]
    GROUND = "vss"

    # ---------- (1) TOPOLOGY ----------
    MOSFETS = [
        Instance("M1", "nmos", d="IREF", g="IREF", s="gnd", b="gnd"),
        Instance("M2", "nmos", d="VOUT", g="IREF", s="gnd", b="gnd"),
    ]

    PASSIVES = [
        Passive("CL", "cap", a="VOUT", b="gnd", external=True),
    ]

    VSOURCES = [
        VSource("VDD", p="vdd", n="gnd", supply=True),
    ]

    SIGNAL_NODES = set()

    # ---------- (2) KNOBS ----------
    KNOBS = [
        Knob("M1_GMID", role="op", sets_width_of="M1"),
        Knob("M1_L", role="geom"),
    ]

    RECORNER_RESOLVE = []

    # ---------- (3) UNKNOWNS ----------
    UNKNOWNS = [
        Unknown("M1_VDS", seed=lambda c: c["vdd"] / 3, bound=lambda c: (0.02, c["vdd"])),
    ]

    # ---------- (4) SOLVE_POINT ----------
    def solve_point(self, v, dev, cond):
        VOUT_DC = cond["vout_dc"]
        IREF = cond["iref"]
        VDD = cond["vdd"]
        K = cond["k"]

        M1 = dev.nmos(gmid=v.M1_GMID, L=v.M1_L, vds=v.M1_VDS, vsb=0.0)
        M2 = dev.nmos(gmid=v.M1_GMID, L=v.M1_L, vds=VOUT_DC, vsb=0.0)

        ID = {
            "M1": IREF,
            "M2": K * IREF,
        }

        W = {
            "M1": ID["M1"] / M1.jd,
            "M2": ID["M2"] / M2.jd,
        }

        L = {"M1": v.M1_L, "M2": v.M1_L}
        GMID = {"M1": v.M1_GMID, "M2": v.M1_GMID}

        ntw = dev.nmos.table_width
        ss = {
            "M1": M1.small_signal(v.M1_GMID, ID["M1"], W["M1"], ntw, use_gmb=False),
            "M2": M2.small_signal(v.M1_GMID, ID["M2"], W["M2"], ntw, use_gmb=False),
        }

        return State(
            M1=M1,
            M2=M2,
            W=W,
            L=L,
            ID=ID,
            GMID=GMID,
            ss=ss,
            VDD=VDD,
            IREF=IREF,
            K=K,
            VOUT_DC=VOUT_DC,
            M1_VDS=v.M1_VDS,
        )

    # ---------- (5) RESIDUALS ----------
    def residuals(self, b):
        return [
            vres(b.M1_VDS, b.M1.vgs),
        ]

    # ---------- (7) SPECS ----------
    def specs(self, b, cond):
        Rout = 1.0 / (b.M2.gds_id * b.ID["M2"])
        Area = b.L["M1"] * b.W["M1"] + b.L["M2"] * b.W["M2"]
        Itotal = b.ID["M1"] + b.ID["M2"]
        Vcompliance = b.M2.vdsat

        return {
            "Rout": Rout,
            "Vcompliance": Vcompliance,
            "Iout": b.ID["M2"],
            "Area": Area,
            "Itotal": Itotal,
        }

    # ---------- (8) NETLIST HOOKS ----------
    def netlist_context(self, corner, ref_op=None) -> dict:
        return {"iref": corner.cond("iref"), "k": corner.cond("k")}


__all__ = ["Circuit", "run", "Knob"]
