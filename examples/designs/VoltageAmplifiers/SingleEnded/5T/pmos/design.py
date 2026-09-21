from __future__ import annotations

import numpy as np

from mosplot.optimizer import (
    CircuitModel,
    Instance,
    Knob,
    Passive,
    Spec,
    State,
    Unknown,
    VSource,
    build_ss_model,
    run,
    vres,
    rres,
)


class Circuit(CircuitModel):
    NAME = "amp"
    PORTS = ["VINN", "VINP", "VOUT", "vdd", "vss"]
    GROUND = "vss"

    # ---------- (1) TOPOLOGY ----------
    MOSFETS = [
        Instance("M1a", "pmos", d="n01", g="VINP", s="n02", b="vdd"),
        Instance("M1b", "pmos", d="VOUT", g="VINN", s="n02", b="vdd"),
        Instance("M2a", "nmos", d="n01", g="n01", s="vss", b="vss"),
        Instance("M2b", "nmos", d="VOUT", g="n01", s="vss", b="vss"),
        Instance("M3", "pmos", d="n02", g="vbp", s="vdd", b="vdd"),
    ]
    PASSIVES = [
        Passive("COUT", "cap", a="VOUT", b="vss", external=True),
    ]
    VSOURCES = [
        VSource("VDD", p="vdd", n="vss", supply=True),
        VSource("VBP", p="vbp", n="vss", mirror="M3"),
    ]
    SIGNAL_NODES = {"VINP", "VINN"}

    # ---------- (2) KNOBS (names + roles only; bounds live in the config) ----------
    KNOBS = [
        Knob("M1a_GMID", role="op", sets_width_of="M1a"),
        Knob("M2a_GMID", role="op", sets_width_of="M2a"),
        Knob("M3_GMID", role="op", sets_width_of="M3"),
        Knob("M1a_L", role="geom"),
        Knob("M2a_L", role="geom"),
        Knob("M3_L", role="geom"),
        Knob("M1a_ID", role="external"),
    ]
    # Multicorner [1]: the diode reference currents IREF are CONSERVED across corners
    # (freeze_extra + recorner_residuals below); the branch currents are re-solved so the
    # tail's delivered current drifts as the real fixed-IREF mirror makes it.
    RECORNER_RESOLVE = ["M1a_ID"]

    # ---------- (3) UNKNOWNS (self-referential lookup inputs) ----------
    UNKNOWNS = [
        Unknown("M1a_VDS", seed=lambda c: c["vdd"] / 3, bound=lambda c: (0.02, c["vdd"])),
        Unknown("M2a_VDS", seed=lambda c: c["vdd"] / 3, bound=lambda c: (0.02, c["vdd"])),
        Unknown("Mvbp_GMID", seed=lambda c: 12.0, bound=lambda c: (5.0, 25.0)),  # diode replica gm/ID [1]
    ]

    # ---------- (4) SOLVE_POINT ----------
    def solve_point(self, v: State, dev, cond) -> State:
        COUT = cond["cout"]
        VDD = cond["vdd"]
        VIN_CM = cond["vin_cm"]
        VOUT_DC = cond["vout_dc"]

        M2a = dev.nmos(gmid=v.M2a_GMID, L=v.M2a_L, vds=v.M2a_VDS, vsb=0.0)
        M1a = dev.pmos(gmid=v.M1a_GMID, L=v.M1a_L, vds=v.M1a_VDS, vsb=0.0)
        M1a_VDS = v.M1a_VDS
        M1b_VDS = M1a.vds_used
        M1b = dev.pmos(gmid=v.M1a_GMID, L=v.M1a_L, vds=M1b_VDS, vsb=0.0)
        M3_VDS = VDD - VIN_CM - M1a.vgs
        M3 = dev.pmos(gmid=v.M3_GMID, L=v.M3_L, vds=M3_VDS, vsb=0.0)
        M2b_VDS = VIN_CM + M1a.vgs - M1b.vds_used
        M2b = dev.nmos(gmid=v.M2a_GMID, L=v.M2a_L, vds=M2b_VDS, vsb=0.0)
        Mvbp = dev.pmos(gmid=v.Mvbp_GMID, L=v.M3_L, vds=M3.vgs, vsb=0.0)

        ID = {}
        ID["M1a"] = v.M1a_ID
        ID["M1b"] = ID["M1a"]
        ID["M2a"] = ID["M1a"]
        ID["M2b"] = ID["M2a"]
        ID["M3"] = ID["M1a"] + ID["M1b"]
        ID["Mvbp"] = ID["M3"]
        IREF_Mvbp = ID["M3"] * Mvbp.jd / M3.jd

        W = {}
        W["M1a"] = ID["M1a"] / M1a.jd
        W["M1b"] = W["M1a"]
        W["M2a"] = ID["M2a"] / M2a.jd
        W["M2b"] = W["M2a"]
        W["M3"] = ID["M3"] / M3.jd
        W["Mvbp"] = W["M3"]  # 1:1 bias diode replica

        L = {
            "M1a": v.M1a_L,
            "M1b": v.M1a_L,
            "M2a": v.M2a_L,
            "M2b": v.M2a_L,
            "M3": v.M3_L,
            "Mvbp": v.M3_L,
        }
        GMID = {
            "M1a": v.M1a_GMID,
            "M1b": v.M1a_GMID,
            "M2a": v.M2a_GMID,
            "M2b": v.M2a_GMID,
            "M3": v.M3_GMID,
            "Mvbp": v.M3_GMID,
        }

        ptw, ntw = dev.pmos.table_width, dev.nmos.table_width
        ss = {
            "M1a": M1a.small_signal(v.M1a_GMID, ID["M1a"], W["M1a"], ptw, use_gmb=False),
            "M1b": M1b.small_signal(v.M1a_GMID, ID["M1b"], W["M1b"], ptw, use_gmb=False),
            "M2a": M2a.small_signal(v.M2a_GMID, ID["M2a"], W["M2a"], ntw, use_gmb=False),
            "M2b": M2b.small_signal(v.M2a_GMID, ID["M2b"], W["M2b"], ntw, use_gmb=False),
            "M3": M3.small_signal(v.M3_GMID, ID["M3"], W["M3"], ptw, use_gmb=False),
        }

        VBP = VDD - M3.vgs

        return State(
            W=W,
            L=L,
            ID=ID,
            GMID=GMID,
            ss=ss,
            M1a=M1a,
            M1b=M1b,
            M2a=M2a,
            M3=M3,
            M2b=M2b,
            Mvbp=Mvbp,
            IREF_Mvbp=IREF_Mvbp,
            COUT=COUT,
            VDD=VDD,
            VIN_CM=VIN_CM,
            VBP=VBP,
            VOUT_DC=VOUT_DC,
            M1a_VDS=M1a_VDS,
            M1b_VDS=M1b_VDS,
            M2a_VDS=v.M2a_VDS,
            M3_VDS=M3_VDS,
            M2b_VDS=M2b_VDS,
            Mvbp_GMID=v.Mvbp_GMID,
        )

    # ---------- (5) RESIDUALS (node closures the solver drives to zero) ----------
    def residuals(self, b) -> list:
        return [
            vres(b.M1a_VDS, b.VIN_CM + b.M1a.vgs - b.M2a.vgs),
            vres(b.M2a_VDS, b.M2a.vgs),
            vres(b.Mvbp.vgs, b.M3.vgs),  # diode replica shares the master's VGS
        ]

    # ---------- (7) MULTICORNER CONSERVATION ----------
    def freeze_extra(self, b) -> dict:
        return {"IREF_Mvbp": b.IREF_Mvbp}

    def recorner_residuals(self, b, frozen) -> list:
        e = frozen["extra"]
        return [rres(b.IREF_Mvbp, e["IREF_Mvbp"], e["IREF_Mvbp"])]

    # ---------- (6) SPECS ----------
    def specs(self, b, cond) -> dict:
        ss = build_ss_model(
            self.MOSFETS,
            self.PASSIVES,
            self.VSOURCES,
            b.ss,
            {"COUT": cond["cout"]},
            signal_nodes=self.SIGNAL_NODES,
        )
        ac = ss.transfer(inputs={"VINP": 0.5, "VINN": -0.5}, output={"VOUT": 1.0})
        ac_cm = ss.transfer(inputs={"VINP": 1.0, "VINN": 1.0}, output={"VOUT": 1.0})
        out = {
            "GBW": ac.ugf(),
            "AC Gain (dB)": 20.0 * np.log10(max(abs(ac.gain()), 1e-300)),
            "PM": ac.phase_margin(),
            "DC CMR (dB)": 20.0 * np.log10(max(ac.rejection(ac_cm), 1e-300)),
        }
        Area = (
            b.L["M1a"] * b.W["M1a"]
            + b.L["M1b"] * b.W["M1b"]
            + b.L["M2a"] * b.W["M2a"]
            + b.L["M2b"] * b.W["M2b"]
            + b.L["M3"] * b.W["M3"]
            + b.L["Mvbp"] * b.W["Mvbp"]
        )
        out["Area"] = Area
        Itotal = b.ID["M1a"] + b.ID["M1b"]
        out["Itotal"] = Itotal
        VBP = b.VDD - b.M3.vgs
        out["VBP"] = VBP
        VOUT_DC = b.VIN_CM + b.M1a.vgs - b.M1b.vds_used
        out["VOUT_DC"] = VOUT_DC
        VOUT_MAX = -b.M1b.vdsat + b.M1b.vgs + b.VIN_CM
        out["VOUT_MAX"] = VOUT_MAX
        VOUT_MIN = b.M2b.vdsat
        out["VOUT_MIN"] = VOUT_MIN
        VIN_MAX = -b.M1a.vgs - b.M3.vdsat + b.VDD
        out["VIN_MAX"] = VIN_MAX
        VIN_MIN = max(
            b.M1a.vdsat - b.M1b.vgs + b.M2b.vgs,
            b.M1b.vdsat - b.M1b.vgs + VOUT_DC,
        )
        out["VIN_MIN"] = VIN_MIN
        Output_Swing = VOUT_MAX - VOUT_MIN
        out["Output_Swing"] = Output_Swing
        return out

    # ---------- netlist hooks ----------
    def mirror_currents(self, ref_op) -> dict:
        return {"VBP": ref_op.IREF_Mvbp}

    def netlist_context(self, corner, ref_op=None) -> dict:
        return {"vcm": corner.cond("vin_cm"), "vout_dc": corner.cond("vout_dc")}


__all__ = ["Circuit", "run", "Knob", "Spec"]
