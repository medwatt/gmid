from __future__ import annotations

from mosplot.optimizer import (
    CircuitModel,
    Instance,
    Knob,
    Passive,
    State,
    Unknown,
    VSource,
    build_ss_model,
    run,
    vres,
)


class Circuit(CircuitModel):
    NAME = "cascode_mirror"
    PORTS = ["IREF", "VOUT", "vdd", "vss"]
    GROUND = "vss"

    MOSFETS = [
        Instance("M1", "pmos", d="n01", g="n01", s="vdd", b="vdd"),
        Instance("M4", "pmos", d="IREF", g="IREF", s="n01", b="vdd"),
        Instance("M2", "pmos", d="n02", g="n01", s="vdd", b="vdd"),
        Instance("M3", "pmos", d="VOUT", g="IREF", s="n02", b="vdd"),
    ]

    PASSIVES = [
        Passive("CL", "cap", a="VOUT", b="gnd", external=True),
    ]

    VSOURCES = [
        VSource("VDD", p="vdd", n="gnd", supply=True),
    ]

    SIGNAL_NODES = set()

    KNOBS = [
        Knob("M1_GMID", role="op", sets_width_of="M1"),
        Knob("M1_L", role="geom"),
        Knob("M3_GMID", role="op", sets_width_of="M4"),
        Knob("M3_L", role="geom"),
    ]

    RECORNER_RESOLVE = []

    UNKNOWNS = [
        Unknown("M1_VDS", seed=lambda c: c["vdd"] / 4, bound=lambda c: (0.02, c["vdd"])),
        Unknown("M4_VDS", seed=lambda c: c["vdd"] / 4, bound=lambda c: (0.02, c["vdd"])),
        Unknown("M2_VDS", seed=lambda c: c["vdd"] / 4, bound=lambda c: (0.02, c["vdd"])),
    ]

    def solve_point(self, v, dev, cond):
        VOUT_DC = cond["vout_dc"]
        IREF = cond["iref"]
        VDD = cond["vdd"]
        K = cond["k"]

        M1 = dev.pmos(gmid=v.M1_GMID, L=v.M1_L, vds=v.M1_VDS, vsb=0.0)
        M4 = dev.pmos(gmid=v.M3_GMID, L=v.M3_L, vds=v.M4_VDS, vsb=v.M1_VDS)

        M3_VDS = VDD - v.M2_VDS - VOUT_DC
        M3 = dev.pmos(gmid=v.M3_GMID, L=v.M3_L, vds=M3_VDS, vsb=v.M2_VDS)
        M2 = dev.pmos(gmid=v.M1_GMID, L=v.M1_L, vds=v.M2_VDS, vsb=0.0)

        ID = {
            "M1": IREF,
            "M4": IREF,
            "M2": K * IREF,
            "M3": K * IREF,
        }

        W = {
            "M1": ID["M1"] / M1.jd,
            "M4": ID["M4"] / M4.jd,
            "M2": ID["M2"] / M2.jd,
            "M3": ID["M3"] / M3.jd,
        }

        L = {
            "M1": v.M1_L,
            "M4": v.M3_L,
            "M2": v.M1_L,
            "M3": v.M3_L,
        }

        GMID = {
            "M1": v.M1_GMID,
            "M4": v.M3_GMID,
            "M2": v.M1_GMID,
            "M3": v.M3_GMID,
        }

        ptw = dev.pmos.table_width
        ss = {
            "M1": M1.small_signal(v.M1_GMID, ID["M1"], W["M1"], ptw),
            "M4": M4.small_signal(v.M3_GMID, ID["M4"], W["M4"], ptw),
            "M2": M2.small_signal(v.M1_GMID, ID["M2"], W["M2"], ptw),
            "M3": M3.small_signal(v.M3_GMID, ID["M3"], W["M3"], ptw),
        }

        return State(
            M1=M1,
            M4=M4,
            M2=M2,
            M3=M3,
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
            M4_VDS=v.M4_VDS,
            M2_VDS=v.M2_VDS,
        )

    def residuals(self, b):
        return [
            vres(b.M1_VDS, b.M1.vgs),
            vres(b.M4_VDS, b.M4.vgs),
            vres(b.M2_VDS, b.M1.vgs + b.M4.vgs - b.M3.vgs),
        ]

    def specs(self, b, cond):
        ss_model = build_ss_model(
            self.MOSFETS,
            self.PASSIVES,
            self.VSOURCES,
            b.ss,
            {"CL": cond.get("cout", 1e-12)},
        )
        Rout = ss_model.port("VOUT").resistance()

        Area = (
            b.L["M1"] * b.W["M1"]
            + b.L["M2"] * b.W["M2"]
            + b.L["M3"] * b.W["M3"]
            + b.L["M4"] * b.W["M4"]
        )
        Itotal = b.ID["M1"] + b.ID["M2"]
        Vcompliance = b.M2_VDS + b.M3.vdsat

        return {
            "Rout": Rout,
            "Vcompliance": Vcompliance,
            "Iout": b.ID["M2"],
            "Area": Area,
            "Itotal": Itotal,
        }

    def netlist_context(self, corner, ref_op=None) -> dict:
        return {"iref": corner.cond("iref"), "k": corner.cond("k")}


__all__ = ["Circuit", "run", "Knob"]
