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
        Instance("M1a", "nmos", d="n01", g="VINN", s="n02", b="vss"),
        Instance("M1b", "nmos", d="n03", g="VINP", s="n02", b="vss"),
        Instance("M2a", "pmos", d="n01", g="n01", s="vdd", b="vdd"),
        Instance("M2b", "pmos", d="n03", g="n01", s="vdd", b="vdd"),
        Instance("M3", "nmos", d="n02", g="vbn", s="vss", b="vss"),
        Instance("M4", "pmos", d="VOUT", g="n03", s="vdd", b="vdd"),
        Instance("M5", "nmos", d="VOUT", g="vbn", s="vss", b="vss"),
    ]
    PASSIVES = [
        Passive("Rz", "res", a="n04", b="n03"),
        Passive("CC", "cap", a="n04", b="VOUT"),
        Passive("COUT", "cap", a="VOUT", b="vss", external=True),
    ]
    VSOURCES = [
        VSource("VDD", p="vdd", n="vss", supply=True),
        VSource("VBN", p="vbn", n="vss", mirror="M3"),
    ]
    SIGNAL_NODES = {"VINP", "VINN"}

    # ---------- (2) KNOBS ----------
    KNOBS = [
        Knob("M1a_GMID", role="op", sets_width_of="M1a"),
        Knob("M2a_GMID", role="op", sets_width_of="M2a"),
        Knob("M3_GMID", role="op", sets_width_of="M3"),
        Knob("M4_GMID", role="op", sets_width_of="M4"),
        Knob("M1a_L", role="geom"),
        Knob("M2a_L", role="geom"),
        Knob("M3_L", role="geom"),
        Knob("M4_L", role="geom"),
        Knob("M5_over_M3", role="geom"),
        Knob("CC", role="geom"),
        Knob("Rz", role="geom"),
        Knob("M1a_ID", role="external"),
    ]
    RECORNER_RESOLVE = ["M1a_ID"]

    # ---------- (3) UNKNOWNS ----------
    UNKNOWNS = [
        Unknown("M1b_VDS", seed=lambda c: c["vdd"] / 3, bound=lambda c: (0.02, c["vdd"])),
        Unknown("M1b_VSB", seed=lambda c: 0, bound=lambda c: (0.02, c["vdd"])),
        Unknown("M2a_VDS", seed=lambda c: c["vdd"] / 3, bound=lambda c: (0.02, c["vdd"])),
        Unknown("Mvbn_GMID", seed=lambda c: 12.0, bound=lambda c: (4.0, 30.0)),  # diode replica gm/ID [1]
    ]

    # ---------- (4) SOLVE_POINT ----------
    def solve_point(self, v: State, dev, cond) -> State:
        COUT = cond["cout"]
        VDD = cond["vdd"]
        VIN_CM = cond["vin_cm"]
        VOUT_DC = cond["vout_dc"]

        M4_VDS = VDD - VOUT_DC
        M4 = dev.pmos(gmid=v.M4_GMID, L=v.M4_L, vds=M4_VDS, vsb=0.0)
        M2a = dev.pmos(gmid=v.M2a_GMID, L=v.M2a_L, vds=v.M2a_VDS, vsb=0.0)
        M2b_VDS = M4.vgs
        M2b = dev.pmos(gmid=v.M2a_GMID, L=v.M2a_L, vds=M2b_VDS, vsb=0.0)
        M1b = dev.nmos(gmid=v.M1a_GMID, L=v.M1a_L, vds=v.M1b_VDS, vsb=v.M1b_VSB)
        M1a_VDS = VDD - VIN_CM + M1b.vgs - M2b.vgs
        M1a_VSB = v.M1b_VSB
        M1a = dev.nmos(gmid=v.M1a_GMID, L=v.M1a_L, vds=M1a_VDS, vsb=M1a_VSB)
        M3_VDS = VIN_CM - M1b.vgs
        M3 = dev.nmos(gmid=v.M3_GMID, L=v.M3_L, vds=M3_VDS, vsb=0.0)
        M5_VDS = VOUT_DC
        M5 = dev.nmos(gmid=v.M3_GMID, L=v.M3_L, vds=M5_VDS, vsb=0.0)
        Mvbn = dev.nmos(gmid=v.Mvbn_GMID, L=v.M3_L, vds=M3.vgs, vsb=0.0)

        ID = {}
        ID["M1a"] = v.M1a_ID
        ID["M1b"] = ID["M1a"]
        ID["M2a"] = ID["M1a"]
        ID["M2b"] = ID["M2a"]
        ID["M3"] = ID["M1a"] + ID["M1b"]
        ID["M5"] = ID["M3"] * v.M5_over_M3
        ID["M4"] = ID["M5"]

        W = {}
        W["M1a"] = ID["M1a"] / M1a.jd
        W["M1b"] = W["M1a"]
        W["M2a"] = ID["M2a"] / M2a.jd
        W["M3"] = ID["M3"] / M3.jd
        W["M4"] = ID["M4"] / M4.jd
        W["M2b"] = W["M2a"]
        W["M5"] = W["M3"] * v.M5_over_M3
        ID["Mvbn"] = ID["M3"]
        W["Mvbn"] = W["M3"]  # 1:1 bias diode replica
        IREF_Mvbn = ID["M3"] * Mvbn.jd / M3.jd

        L = {
            "M1a": v.M1a_L,
            "M1b": v.M1a_L,
            "M2a": v.M2a_L,
            "M2b": v.M2a_L,
            "M3": v.M3_L,
            "M4": v.M4_L,
            "M5": v.M3_L,
            "Mvbn": v.M3_L,
        }
        GMID = {
            "M1a": v.M1a_GMID,
            "M1b": v.M1a_GMID,
            "M2a": v.M2a_GMID,
            "M3": v.M3_GMID,
            "M2b": v.M2a_GMID,
            "M4": v.M4_GMID,
            "M5": v.M3_GMID,
            "Mvbn": v.Mvbn_GMID,
        }

        ptw, ntw = dev.pmos.table_width, dev.nmos.table_width
        ss = {
            "M1a": M1a.small_signal(v.M1a_GMID, ID["M1a"], W["M1a"], ntw, use_gmb=False),
            "M1b": M1b.small_signal(v.M1a_GMID, ID["M1b"], W["M1b"], ntw, use_gmb=False),
            "M2a": M2a.small_signal(v.M2a_GMID, ID["M2a"], W["M2a"], ptw, use_gmb=False),
            "M2b": M2b.small_signal(v.M2a_GMID, ID["M2b"], W["M2b"], ptw, use_gmb=False),
            "M3": M3.small_signal(v.M3_GMID, ID["M3"], W["M3"], ntw, use_gmb=False),
            "M4": M4.small_signal(v.M4_GMID, ID["M4"], W["M4"], ptw, use_gmb=False),
            "M5": M5.small_signal(v.M3_GMID, ID["M5"], W["M5"], ntw, use_gmb=False),
        }

        VBN = M3.vgs

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
            M4=M4,
            M2b=M2b,
            M5=M5,
            Mvbn=Mvbn,
            IREF_Mvbn=IREF_Mvbn,
            COUT=COUT,
            VDD=VDD,
            VIN_CM=VIN_CM,
            VOUT_DC=VOUT_DC,
            VBN=VBN,
            M5_over_M3=v.M5_over_M3,
            CC=v.CC,
            Rz=v.Rz,
            M1a_VDS=M1a_VDS,
            M1a_VSB=M1a_VSB,
            M1b_VDS=v.M1b_VDS,
            M1b_VSB=v.M1b_VSB,
            M2a_VDS=v.M2a_VDS,
            M3_VDS=M3_VDS,
            M4_VDS=M4_VDS,
            M2b_VDS=M2b_VDS,
            M5_VDS=M5_VDS,
            Mvbn_GMID=v.Mvbn_GMID,
        )

    # ---------- (5) RESIDUALS ----------
    def residuals(self, b) -> list:
        return [
            vres(b.M1b_VDS, b.VDD - b.VIN_CM + b.M1b.vgs - b.M4.vgs),
            vres(b.M1b_VSB, b.VIN_CM - b.M1b.vgs),
            vres(b.M2a_VDS, b.M2a.vgs),
            vres(b.Mvbn.vgs, b.M3.vgs),  # diode replica shares the master's VGS
        ]

    # ---------- (6) MULTICORNER CONSERVATION ----------
    def freeze_extra(self, b) -> dict:
        return {"IREF_Mvbn": b.IREF_Mvbn}

    def recorner_residuals(self, b, frozen) -> list:
        e = frozen["extra"]
        return [rres(b.IREF_Mvbn, e["IREF_Mvbn"], e["IREF_Mvbn"])]

    # ---------- (7) SPECS ----------
    def specs(self, b, cond) -> dict:
        ss = build_ss_model(
            self.MOSFETS,
            self.PASSIVES,
            self.VSOURCES,
            b.ss,
            {"Rz": b.Rz, "CC": b.CC, "COUT": cond["cout"]},
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
        out["VOUT_Error"] = abs(b.M4.vgs - b.M2a.vgs)
        Area = (
            b.L["M1a"] * b.W["M1a"]
            + b.L["M1b"] * b.W["M1b"]
            + b.L["M2a"] * b.W["M2a"]
            + b.L["M2b"] * b.W["M2b"]
            + b.L["M3"] * b.W["M3"]
            + b.L["M4"] * b.W["M4"]
            + b.L["M5"] * b.W["M5"]
            + b.L["Mvbn"] * b.W["Mvbn"]
        )
        out["Area"] = Area
        Itotal = b.ID["M1a"] + b.ID["M2b"] + b.ID["M5"]
        out["Itotal"] = Itotal
        VBN = b.M3.vgs
        out["VBN"] = VBN
        VOUT_MAX = -b.M4.vdsat + b.VDD
        out["VOUT_MAX"] = VOUT_MAX
        VOUT_MIN = b.M5.vdsat
        out["VOUT_MIN"] = VOUT_MIN
        VIN_MAX = min(
            -b.M1a.vdsat + b.M1b.vgs - b.M2b.vgs + b.VDD,
            -b.M1b.vdsat + b.M1b.vgs - b.M4.vgs + b.VDD,
        )
        out["VIN_MAX"] = VIN_MAX
        VIN_MIN = b.M1b.vgs + b.M3.vdsat
        out["VIN_MIN"] = VIN_MIN
        Output_Swing = VOUT_MAX - VOUT_MIN
        out["Output_Swing"] = Output_Swing
        return out

    # ---------- (8) NETLIST HOOKS ----------
    def mirror_currents(self, ref_op) -> dict:
        return {"VBN": ref_op.IREF_Mvbn}

    def passive_values(self, ref_op) -> dict:
        return {"Rz": ref_op.Rz, "CC": ref_op.CC}

    def netlist_context(self, corner, ref_op=None) -> dict:
        return {"vcm": corner.cond("vin_cm"), "vout_dc": corner.cond("vout_dc")}


__all__ = ["Circuit", "run", "Knob", "Spec"]
