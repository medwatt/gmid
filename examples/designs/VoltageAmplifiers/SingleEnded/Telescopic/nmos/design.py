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
        Instance("M1a", "nmos", d="n01", g="VINP", s="n02", b="vss"),
        Instance("M1b", "nmos", d="n03", g="VINN", s="n02", b="vss"),
        Instance("M2a", "nmos", d="n04", g="vcascn", s="n01", b="vss"),
        Instance("M2b", "nmos", d="VOUT", g="vcascn", s="n03", b="vss"),
        Instance("M3a", "pmos", d="n04", g="vcascp", s="n05", b="vdd"),
        Instance("M3b", "pmos", d="VOUT", g="vcascp", s="n06", b="vdd"),
        Instance("M4a", "pmos", d="n05", g="n04", s="vdd", b="vdd"),
        Instance("M4b", "pmos", d="n06", g="n04", s="vdd", b="vdd"),
        Instance("M5", "nmos", d="n02", g="vbn", s="vss", b="vss"),
    ]
    PASSIVES = [
        Passive("COUT", "cap", a="VOUT", b="vss", external=True),
    ]
    VSOURCES = [
        VSource("VDD", p="vdd", n="vss", supply=True),
        VSource("VCASCN", p="vcascn", n="vss"),
        VSource("VCASCP", p="vcascp", n="vss"),
        VSource("VBN", p="vbn", n="vss", mirror="M5"),
    ]
    SIGNAL_NODES = {"VINP", "VINN"}

    # ---------- (2) KNOBS (names + roles only; bounds live in the config) ----------
    KNOBS = [
        Knob("M1a_GMID", role="op", sets_width_of="M1a"),
        Knob("M2a_GMID", role="op", sets_width_of="M2a"),
        Knob("M3a_GMID", role="op", sets_width_of="M3a"),
        Knob("M4a_GMID", role="op", sets_width_of="M4a"),
        Knob("M5_GMID", role="op", sets_width_of="M5"),
        Knob("M1a_L", role="geom"),
        Knob("M2a_L", role="geom"),
        Knob("M3a_L", role="geom"),
        Knob("M4a_L", role="geom"),
        Knob("M5_L", role="geom"),
        Knob("M1a_ID", role="external"),
        Knob("M4a_VDSAT_MARGIN", role="external"),
        Knob("M1a_VDSAT_MARGIN", role="external"),
    ]
    # Multicorner [1]: what stays fixed on the chip is CONSERVED across corners (freeze_extra +
    # recorner_residuals below): the diode reference current, which re-solves the branch
    # current M1a_ID, and the two cascode bias voltages (fixed sources in the netlist), which
    # re-solve the saturation margins they are derived from (M1a for VCASC at the input side,
    # M4a for the load side).
    RECORNER_RESOLVE = ["M1a_ID", "M1a_VDSAT_MARGIN", "M4a_VDSAT_MARGIN"]

    # ---------- (3) UNKNOWNS (node VDS/VSB the DC solver finds) ----------
    # REVIEW: every device VDS is an unknown for robustness. Many are explicit node
    # closures (RHS has no dependence on the device's own VGS) and can be folded into
    # solve_point by hand to shrink the system; the diode/self-referential ones must stay.
    UNKNOWNS = [
        Unknown("M1a_VDS", seed=lambda c: c["vdd"] / 5, bound=lambda c: (0.02, c["vdd"])),
        Unknown("M1a_VSB", seed=lambda c: 0, bound=lambda c: (0.02, c["vdd"])),
        Unknown("M2a_VDS", seed=lambda c: c["vdd"] / 5, bound=lambda c: (0.02, c["vdd"])),
        Unknown("M3a_VDS", seed=lambda c: c["vdd"] / 5, bound=lambda c: (0.02, c["vdd"])),
        Unknown("M4a_VDS", seed=lambda c: c["vdd"] / 5, bound=lambda c: (0.02, c["vdd"])),
        Unknown("M2b_VSB", seed=lambda c: 0, bound=lambda c: (0.02, c["vdd"])),
        Unknown("M3b_VDS", seed=lambda c: c["vdd"] / 5, bound=lambda c: (0.02, c["vdd"])),
        Unknown("M3b_VSB", seed=lambda c: 0, bound=lambda c: (0.02, c["vdd"])),
        Unknown("Mvbn_GMID", seed=lambda c: 12.0, bound=lambda c: (4.0, 30.0)),  # diode replica gm/ID [1]
    ]

    # ---------- (4) SOLVE_POINT ----------
    def solve_point(self, v: State, dev, cond) -> State:
        COUT = cond["cout"]
        VDD = cond["vdd"]
        VIN_CM = cond["vin_cm"]

        M1a = dev.nmos(gmid=v.M1a_GMID, L=v.M1a_L, vds=v.M1a_VDS, vsb=v.M1a_VSB)
        M4a = dev.pmos(gmid=v.M4a_GMID, L=v.M4a_L, vds=v.M4a_VDS, vsb=0.0)
        M5_VDS = VIN_CM - M1a.vgs
        M5 = dev.nmos(gmid=v.M5_GMID, L=v.M5_L, vds=M5_VDS, vsb=0.0)
        M2a_VSB = VIN_CM + M1a.vds_used - M1a.vgs
        M2a = dev.nmos(gmid=v.M2a_GMID, L=v.M2a_L, vds=v.M2a_VDS, vsb=M2a_VSB)
        M3a_VSB = M4a.vds_used
        M3a = dev.pmos(gmid=v.M3a_GMID, L=v.M3a_L, vds=v.M3a_VDS, vsb=M3a_VSB)
        M2b_VDS = M2a.vds_used
        M2b = dev.nmos(gmid=v.M2a_GMID, L=v.M2a_L, vds=M2b_VDS, vsb=v.M2b_VSB)
        M3b = dev.pmos(gmid=v.M3a_GMID, L=v.M3a_L, vds=v.M3b_VDS, vsb=v.M3b_VSB)
        M4b_VDS = M3a.vgs - M3b.vgs + M4a.vds_used
        M4b = dev.pmos(gmid=v.M4a_GMID, L=v.M4a_L, vds=M4b_VDS, vsb=0.0)
        M1b_VDS = M1a.vds_used + M2a.vgs - M2b.vgs
        M1b_VSB = VIN_CM - M1a.vgs
        M1b = dev.nmos(gmid=v.M1a_GMID, L=v.M1a_L, vds=M1b_VDS, vsb=M1b_VSB)
        Mvbn = dev.nmos(gmid=v.Mvbn_GMID, L=v.M5_L, vds=M5.vgs, vsb=0.0)

        ID = {}
        ID["M1a"] = v.M1a_ID
        ID["M1b"] = ID["M1a"]
        ID["M2a"] = ID["M1a"]
        ID["M3a"] = ID["M1a"]
        ID["M4a"] = ID["M1a"]
        ID["M5"] = ID["M1a"] + ID["M1b"]
        ID["M2b"] = ID["M2a"]
        ID["M3b"] = ID["M3a"]
        ID["M4b"] = ID["M4a"]

        W = {}
        W["M1a"] = ID["M1a"] / M1a.jd
        W["M1b"] = W["M1a"]
        W["M2a"] = ID["M2a"] / M2a.jd
        W["M3a"] = ID["M3a"] / M3a.jd
        W["M4a"] = ID["M4a"] / M4a.jd
        W["M5"] = ID["M5"] / M5.jd
        W["M2b"] = W["M2a"]
        W["M3b"] = W["M3a"]
        W["M4b"] = W["M4a"]
        W["Mvbn"] = W["M5"]  # 1:1 bias diode replica
        ID["Mvbn"] = ID["M5"]
        IREF_Mvbn = ID["M5"] * Mvbn.jd / M5.jd

        L = {
            "M1a": v.M1a_L,
            "M1b": v.M1a_L,
            "M2a": v.M2a_L,
            "M3a": v.M3a_L,
            "M4a": v.M4a_L,
            "M5": v.M5_L,
            "M2b": v.M2a_L,
            "M3b": v.M3a_L,
            "M4b": v.M4a_L,
            "Mvbn": v.M5_L,
        }
        GMID = {
            "M1a": v.M1a_GMID,
            "M1b": v.M1a_GMID,
            "M2a": v.M2a_GMID,
            "M3a": v.M3a_GMID,
            "M4a": v.M4a_GMID,
            "M5": v.M5_GMID,
            "M2b": v.M2a_GMID,
            "M3b": v.M3a_GMID,
            "M4b": v.M4a_GMID,
            "Mvbn": v.M5_GMID,
        }

        ptw, ntw = dev.pmos.table_width, dev.nmos.table_width
        ss = {
            "M1a": M1a.small_signal(v.M1a_GMID, ID["M1a"], W["M1a"], ntw, use_gmb=False),
            "M1b": M1b.small_signal(v.M1a_GMID, ID["M1b"], W["M1b"], ntw, use_gmb=False),
            "M2a": M2a.small_signal(v.M2a_GMID, ID["M2a"], W["M2a"], ntw, use_gmb=False),
            "M3a": M3a.small_signal(v.M3a_GMID, ID["M3a"], W["M3a"], ptw, use_gmb=False),
            "M4a": M4a.small_signal(v.M4a_GMID, ID["M4a"], W["M4a"], ptw, use_gmb=False),
            "M5": M5.small_signal(v.M5_GMID, ID["M5"], W["M5"], ntw, use_gmb=False),
            "M2b": M2b.small_signal(v.M2a_GMID, ID["M2b"], W["M2b"], ntw, use_gmb=False),
            "M3b": M3b.small_signal(v.M3a_GMID, ID["M3b"], W["M3b"], ptw, use_gmb=False),
            "M4b": M4b.small_signal(v.M4a_GMID, ID["M4b"], W["M4b"], ptw, use_gmb=False),
        }

        VCASCN = VIN_CM + M1a.vds_used - M1a.vgs + M2a.vgs
        VCASCP = VDD - M3a.vgs - M4a.vds_used
        VBN = M5.vgs
        return State(
            W=W,
            L=L,
            ID=ID,
            GMID=GMID,
            ss=ss,
            M1a=M1a,
            M1b=M1b,
            M2a=M2a,
            M3a=M3a,
            M4a=M4a,
            M5=M5,
            M2b=M2b,
            M3b=M3b,
            M4b=M4b,
            Mvbn=Mvbn,
            IREF_Mvbn=IREF_Mvbn,
            COUT=COUT,
            VDD=VDD,
            VIN_CM=VIN_CM,
            VCASCN=VCASCN,
            VCASCP=VCASCP,
            VBN=VBN,
            M4a_VDSAT_MARGIN=v.M4a_VDSAT_MARGIN,
            M1a_VDSAT_MARGIN=v.M1a_VDSAT_MARGIN,
            M1a_VDS=v.M1a_VDS,
            M1a_VSB=v.M1a_VSB,
            M1b_VDS=M1b_VDS,
            M1b_VSB=M1b_VSB,
            M2a_VDS=v.M2a_VDS,
            M2a_VSB=M2a_VSB,
            M3a_VDS=v.M3a_VDS,
            M3a_VSB=M3a_VSB,
            M4a_VDS=v.M4a_VDS,
            M5_VDS=M5_VDS,
            M2b_VDS=M2b_VDS,
            M2b_VSB=v.M2b_VSB,
            M3b_VDS=v.M3b_VDS,
            M3b_VSB=v.M3b_VSB,
            M4b_VDS=M4b_VDS,
            Mvbn_GMID=v.Mvbn_GMID,
        )

    # ---------- (5) RESIDUALS (node closures the solver drives to zero) ----------
    def residuals(self, b) -> list:
        return [
            vres(b.M1a_VDS, b.M1a.vdsat + b.M1a_VDSAT_MARGIN),
            vres(b.M1a_VSB, b.VIN_CM - b.M1a.vgs),
            vres(b.M2a_VDS, b.VDD - b.VIN_CM - b.M1a.vds_used + b.M1a.vgs - b.M4a.vgs),
            vres(b.M3a_VDS, -b.M4a.vds_used + b.M4a.vgs),
            vres(b.M4a_VDS, b.M4a.vdsat + b.M4a_VDSAT_MARGIN),
            vres(b.M2b_VSB, b.VIN_CM + b.M1a.vds_used - b.M1a.vgs + b.M2a.vgs - b.M2b.vgs),
            vres(
                b.M3b_VDS,
                b.VDD
                - b.VIN_CM
                - b.M1a.vds_used
                + b.M1a.vgs
                - b.M2a.vgs
                - b.M2b.vds_used
                + b.M2b.vgs
                - b.M3a.vgs
                + b.M3b.vgs
                - b.M4a.vds_used,
            ),
            vres(b.M3b_VSB, b.M3a.vgs - b.M3b.vgs + b.M4a.vds_used),
            vres(b.Mvbn.vgs, b.M5.vgs),  # diode replica shares the master's VGS
        ]

    # ---------- (7) MULTICORNER CONSERVATION ----------
    def freeze_extra(self, b) -> dict:
        return {"IREF_Mvbn": b.IREF_Mvbn, "VCASCN": b.VCASCN, "VCASCP": b.VCASCP}

    def recorner_residuals(self, b, frozen) -> list:
        e = frozen["extra"]
        return [
            rres(b.IREF_Mvbn, e["IREF_Mvbn"], e["IREF_Mvbn"]),
            vres(b.VCASCN, e["VCASCN"], 0.05),
            vres(b.VCASCP, e["VCASCP"], 0.05),
        ]

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
            + b.L["M3a"] * b.W["M3a"]
            + b.L["M3b"] * b.W["M3b"]
            + b.L["M4a"] * b.W["M4a"]
            + b.L["M4b"] * b.W["M4b"]
            + b.L["M5"] * b.W["M5"]
            + b.L["Mvbn"] * b.W["Mvbn"]
        )
        out["Area"] = Area
        Itotal = b.ID["M1a"] + b.ID["M1b"]
        out["Itotal"] = Itotal
        VCASCN = b.VIN_CM + b.M1a.vds_used - b.M1a.vgs + b.M2a.vgs
        out["VCASCN"] = VCASCN
        VCASCP = b.VDD - b.M3a.vgs - b.M4a.vds_used
        out["VCASCP"] = VCASCP
        VBN = b.M5.vgs
        out["VBN"] = VBN
        VOUT_DC = b.VIN_CM + b.M1a.vds_used - b.M1a.vgs + b.M2a.vgs + b.M2b.vds_used - b.M2b.vgs
        out["VOUT_DC"] = VOUT_DC
        VOUT_MAX = -b.M3b.vdsat + b.M3b.vgs + b.VCASCP
        out["VOUT_MAX"] = VOUT_MAX
        VOUT_MIN = b.M2b.vdsat - b.M2b.vgs + b.VCASCN
        out["VOUT_MIN"] = VOUT_MIN
        VIN_MAX = min(
            -b.M1a.vdsat + b.M1a.vgs - b.M2a.vgs + b.VCASCN,
            b.M1a.vgs - b.M1b.vdsat - b.M2b.vgs + b.VCASCN,
        )
        out["VIN_MAX"] = VIN_MAX
        VIN_MIN = b.M1a.vgs + b.M5.vdsat
        out["VIN_MIN"] = VIN_MIN
        Output_Swing = VOUT_MAX - VOUT_MIN
        out["Output_Swing"] = Output_Swing
        # saturation margins of the devices whose VDS the fixed cascode biases set: re-solved
        # on every non-reference corner, so they show what the corner does to them
        out["M1a margin"] = b.M1a_VDSAT_MARGIN
        out["M4a margin"] = b.M4a_VDSAT_MARGIN
        out["Margin min"] = min(b.M1a_VDSAT_MARGIN, b.M4a_VDSAT_MARGIN)
        return out

    # ---------- netlist hooks ----------
    def vsource_values(self, ref_op, frozen) -> dict:
        return {"VCASCN": ref_op.VCASCN, "VCASCP": ref_op.VCASCP}

    def mirror_currents(self, ref_op) -> dict:
        return {"VBN": ref_op.IREF_Mvbn}

    def netlist_context(self, corner, ref_op=None) -> dict:
        return {"vcm": corner.cond("vin_cm")}


__all__ = ["Circuit", "run", "Knob", "Spec"]
