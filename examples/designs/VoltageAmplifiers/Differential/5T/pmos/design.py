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
    cmfb_loop_sign,
    run,
    vres,
    rres,
)


class Circuit(CircuitModel):
    NAME = "amp"
    PORTS = ["VINN", "VINP", "VOUTP", "VOUTN", "vdd", "vss"]
    GROUND = "vss"

    # ---------- (1) TOPOLOGY ----------
    MOSFETS = [
        Instance("M1a", "pmos", d="VOUTN", g="VINP", s="n01", b="vdd"),
        Instance("M1b", "pmos", d="VOUTP", g="VINN", s="n01", b="vdd"),
        Instance("M2a", "nmos", d="VOUTN", g="vbn", s="vss", b="vss"),
        Instance("M2b", "nmos", d="VOUTP", g="vbn", s="vss", b="vss"),
        Instance("M3", "pmos", d="n01", g="vcmfb", s="vdd", b="vdd"),
    ]
    PASSIVES = [
        Passive("COUTN", "cap", a="VOUTN", b="vss", external=True),
        Passive("COUTP", "cap", a="VOUTP", b="vss", external=True),
    ]
    VSOURCES = [
        VSource("VDD", p="vdd", n="vss", supply=True),
        VSource("VBN", p="vbn", n="vss", mirror="M2a"),
        VSource("VCMFB", p="vcmfb", n="vss", emit=False),
    ]
    SIGNAL_NODES = {"VINP", "VINN"}

    # ---------- (2) KNOBS ----------
    KNOBS = [
        Knob("M1a_GMID", role="op", sets_width_of="M1a"),
        Knob("M2a_GMID", role="op", sets_width_of="M2a"),
        Knob("M3_GMID", role="op", sets_width_of="M3"),
        Knob("M1a_L", role="geom"),
        Knob("M2a_L", role="geom"),
        Knob("M3_L", role="geom"),
        Knob("M1a_ID", role="external"),
        Knob("M3_VDSAT_MARGIN", role="external"),
    ]
    RECORNER_RESOLVE = ["M1a_ID"]

    # ---------- (3) UNKNOWNS ----------
    UNKNOWNS = [
        Unknown("M1a_VDS", seed=lambda c: c["vdd"] / 3, bound=lambda c: (0.02, c["vdd"])),
        Unknown("M2a_VDS", seed=lambda c: c["vdd"] / 3, bound=lambda c: (0.02, c["vdd"])),
        Unknown("M3_VDS", seed=lambda c: c["vdd"] / 3, bound=lambda c: (0.02, c["vdd"])),
        Unknown("Mvbn_GMID", seed=lambda c: 12.0, bound=lambda c: (4.0, 30.0)),    # diode replica gm/ID [1]
        Unknown("Mvcmfb_GMID", seed=lambda c: 12.0, bound=lambda c: (4.0, 30.0)),  # diode replica gm/ID [1]
    ]

    # ---------- (4) SOLVE_POINT ----------
    def solve_point(self, v: State, dev, cond) -> State:
        COUT = cond["cout"]
        VDD = cond["vdd"]
        VIN_CM = cond["vin_cm"]
        VOUT_CM = cond["vout_cm"]

        M1a = dev.pmos(gmid=v.M1a_GMID, L=v.M1a_L, vds=v.M1a_VDS, vsb=0.0)
        M2a = dev.nmos(gmid=v.M2a_GMID, L=v.M2a_L, vds=v.M2a_VDS, vsb=0.0)
        M3 = dev.pmos(gmid=v.M3_GMID, L=v.M3_L, vds=v.M3_VDS, vsb=0.0)
        Mvbn = dev.nmos(gmid=v.Mvbn_GMID, L=v.M2a_L, vds=M2a.vgs, vsb=0.0)
        Mvcmfb = dev.pmos(gmid=v.Mvcmfb_GMID, L=v.M3_L, vds=M3.vgs, vsb=0.0)
        M1b = M1a
        M2b = M2a

        ID = {}
        ID["M1a"] = v.M1a_ID
        ID["M2a"] = ID["M1a"]
        ID["M1b"] = ID["M1a"]
        ID["M2b"] = ID["M2a"]
        ID["M3"] = 2 * ID["M1a"]
        ID["Mvbn"] = ID["M2a"]
        ID["Mvcmfb"] = ID["M3"]
        IREF_Mvbn = ID["M2a"] * Mvbn.jd / M2a.jd
        IREF_Mvcmfb = ID["M3"] * Mvcmfb.jd / M3.jd

        W = {}
        W["M1a"] = ID["M1a"] / M1a.jd
        W["M2a"] = ID["M2a"] / M2a.jd
        W["M1b"] = W["M1a"]
        W["M2b"] = W["M2a"]
        W["M3"] = ID["M3"] / M3.jd
        W["Mvbn"] = W["M2a"]  # 1:1 bias diode replica
        W["Mvcmfb"] = W["M3"]  # 1:1 bias diode replica

        L = {
            "M1a": v.M1a_L,
            "M2a": v.M2a_L,
            "M1b": v.M1a_L,
            "M2b": v.M2a_L,
            "M3": v.M3_L,
            "Mvbn": v.M2a_L,
            "Mvcmfb": v.M3_L,
        }
        GMID = {
            "M1a": v.M1a_GMID,
            "M2a": v.M2a_GMID,
            "M1b": v.M1a_GMID,
            "M2b": v.M2a_GMID,
            "M3": v.M3_GMID,
            "Mvbn": v.M2a_GMID,
            "Mvcmfb": v.M3_GMID,
        }

        ptw, ntw = dev.pmos.table_width, dev.nmos.table_width
        ss = {
            "M1a": M1a.small_signal(v.M1a_GMID, ID["M1a"], W["M1a"], ptw, use_gmb=False),
            "M2a": M2a.small_signal(v.M2a_GMID, ID["M2a"], W["M2a"], ntw, use_gmb=False),
            "M1b": M1b.small_signal(v.M1a_GMID, ID["M1b"], W["M1b"], ptw, use_gmb=False),
            "M2b": M2b.small_signal(v.M2a_GMID, ID["M2b"], W["M2b"], ntw, use_gmb=False),
            "M3": M3.small_signal(v.M3_GMID, ID["M3"], W["M3"], ptw, use_gmb=False),
        }

        VBN = M2a.vgs
        VCMFB = VDD - M3.vgs
        return State(
            W=W,
            L=L,
            ID=ID,
            GMID=GMID,
            ss=ss,
            M1a=M1a,
            M2a=M2a,
            M3=M3,
            M1b=M1b,
            M2b=M2b,
            Mvbn=Mvbn,
            Mvcmfb=Mvcmfb,
            IREF_Mvbn=IREF_Mvbn,
            IREF_Mvcmfb=IREF_Mvcmfb,
            COUT=COUT,
            VDD=VDD,
            VIN_CM=VIN_CM,
            VOUT_CM=VOUT_CM,
            VBN=VBN,
            VCMFB=VCMFB,
            M1a_VDS=v.M1a_VDS,
            M2a_VDS=v.M2a_VDS,
            M3_VDS=v.M3_VDS,
            M3_VDSAT_MARGIN=v.M3_VDSAT_MARGIN,
            Mvbn_GMID=v.Mvbn_GMID,
            Mvcmfb_GMID=v.Mvcmfb_GMID,
        )

    # ---------- (5) RESIDUALS ----------
    def residuals(self, b) -> list:
        return [
            vres(b.M1a_VDS, b.VIN_CM + b.M1a.vgs - b.VOUT_CM),
            vres(b.M2a_VDS, b.VOUT_CM),
            vres(b.M3_VDS, b.M3.vdsat + b.M3_VDSAT_MARGIN),
            vres(b.Mvbn.vgs, b.M2a.vgs),  # diode replica shares the master's VGS
            vres(b.Mvcmfb.vgs, b.M3.vgs), # diode replica shares the master's VGS
        ]

    # ---------- (6) MULTICORNER CONSERVATION ----------
    def freeze_extra(self, b) -> dict:
        return {"IREF_Mvbn": b.IREF_Mvbn}

    def recorner_residuals(self, b, frozen) -> list:
        e = frozen["extra"]
        return [
            rres(b.IREF_Mvbn, e["IREF_Mvbn"], e["IREF_Mvbn"]),
        ]

    # ---------- (7) SPECS ----------
    def specs(self, b, cond) -> dict:
        ss = build_ss_model(
            self.MOSFETS,
            self.PASSIVES,
            self.VSOURCES,
            b.ss,
            {"COUTN": cond["cout"], "COUTP": cond["cout"]},
            signal_nodes=self.SIGNAL_NODES,
        )
        ac = ss.transfer(inputs={"VINP": 0.5, "VINN": -0.5}, output={"VOUTP": 1.0, "VOUTN": -1.0})
        ac_cm = ss.transfer(inputs={"VINP": 1.0, "VINN": 1.0}, output={"VOUTP": 1.0, "VOUTN": -1.0})
        out = {
            "GBW": ac.ugf(),
            "AC Gain (dB)": 20.0 * np.log10(max(abs(ac.gain()), 1e-300)),
            "PM": ac.phase_margin(),
            "DC CMR (dB)": 20.0 * np.log10(max(ac.rejection(ac_cm), 1e-300)),
        }
        Area = (
            2 * b.L["M1a"] * b.W["M1a"]
            + 2 * b.L["M2a"] * b.W["M2a"]
            + b.L["M3"] * b.W["M3"]
            + b.L["Mvbn"] * b.W["Mvbn"]
            + b.L["Mvcmfb"] * b.W["Mvcmfb"]
        )
        out["Area"] = Area
        Itotal = 2 * b.ID["M1a"]
        out["Itotal"] = Itotal
        out["VBN"] = b.VBN
        out["VCMFB"] = b.VCMFB
        VOUT_MAX = b.M1a.vgs - b.M1b.vdsat + b.VIN_CM
        out["VOUT_MAX"] = VOUT_MAX
        VOUT_MIN = b.M2b.vdsat
        out["VOUT_MIN"] = VOUT_MIN
        VIN_MAX = -b.M1a.vgs - b.M3.vdsat + b.VDD
        out["VIN_MAX"] = VIN_MAX
        VIN_MIN = b.M1a.vdsat - b.M1a.vgs + b.VOUT_CM
        out["VIN_MIN"] = VIN_MIN
        Output_Swing = VOUT_MAX - VOUT_MIN
        out["Output_Swing"] = Output_Swing
        return out

    # ---------- (8) NETLIST HOOKS ----------
    def mirror_currents(self, ref_op) -> dict:
        return {"VBN": ref_op.IREF_Mvbn}

    def extra_netlist_lines(self, b) -> list:
        lines = []
        sign = cmfb_loop_sign(self, b, "vcmfb", ["VOUTP", "VOUTN"])
        nominal = float(b.VCMFB)
        lines.append(
            f"Bvcmfb (vcmfb vss) bsource v={nominal:.6g} + ({sign:+.0f})*0.4*tanh(75*(vout_cm - (v(VOUTP)+v(VOUTN))/2))"
        )
        return lines

    def netlist_context(self, corner, ref_op=None) -> dict:
        return {"vcm": corner.cond("vin_cm"), "vout_cm": corner.cond("vout_cm")}


__all__ = ["Circuit", "run", "Knob", "Spec"]
