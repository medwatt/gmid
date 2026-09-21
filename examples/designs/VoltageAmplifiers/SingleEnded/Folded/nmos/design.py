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
        Instance("M2a", "pmos", d="n01", g="vbp", s="vdd", b="vdd"),
        Instance("M2b", "pmos", d="n03", g="vbp", s="vdd", b="vdd"),
        Instance("M3a", "pmos", d="n04", g="vcascp", s="n01", b="vdd"),
        Instance("M3b", "pmos", d="VOUT", g="vcascp", s="n03", b="vdd"),
        Instance("M4a", "nmos", d="n04", g="vcascn", s="n05", b="vss"),
        Instance("M4b", "nmos", d="VOUT", g="vcascn", s="n06", b="vss"),
        Instance("M5a", "nmos", d="n05", g="n04", s="vss", b="vss"),
        Instance("M5b", "nmos", d="n06", g="n04", s="vss", b="vss"),
        Instance("M6", "nmos", d="n02", g="vbn", s="vss", b="vss"),
    ]
    PASSIVES = [
        Passive("COUT", "cap", a="VOUT", b="vss", external=True),
    ]
    VSOURCES = [
        VSource("VDD", p="vdd", n="vss", supply=True),
        VSource("VCASCN", p="vcascn", n="vss"),
        VSource("VCASCP", p="vcascp", n="vss"),
        VSource("VBP", p="vbp", n="vss", mirror="M2a"),
        VSource("VBN", p="vbn", n="vss", mirror="M6"),
    ]
    SIGNAL_NODES = {"VINP", "VINN"}

    # ---------- (2) KNOBS (names + roles only; bounds live in the config) ----------
    KNOBS = [
        Knob("M1a_GMID", role="op", sets_width_of="M1a"),
        Knob("M2a_GMID", role="op", sets_width_of="M2a"),
        Knob("M3a_GMID", role="op", sets_width_of="M3a"),
        Knob("M4a_GMID", role="op", sets_width_of="M4a"),
        Knob("M5a_GMID", role="op", sets_width_of="M5a"),
        Knob("M6_GMID", role="op", sets_width_of="M6"),
        Knob("M1a_L", role="geom"),
        Knob("M2a_L", role="geom"),
        Knob("M3a_L", role="geom"),
        Knob("M4a_L", role="geom"),
        Knob("M5a_L", role="geom"),
        Knob("M6_L", role="geom"),
        Knob("M5a_ID_over_M1a_ID", role="geom"),
        Knob("M1a_ID", role="external"),
        Knob("M5a_VDSAT_MARGIN", role="external"),
        Knob("M2a_VDSAT_MARGIN", role="external"),
    ]
    # Multicorner [1]: what stays fixed on the chip when the process moves is CONSERVED across
    # corners (freeze_extra + recorner_residuals below): the two diode reference currents, and
    # the two cascode bias voltages, which the netlist emits as fixed sources. Each conserved
    # quantity re-solves one knob: the two reference currents fix the tail current (M1a_ID)
    # and the folding current (M5a_ID_over_M1a_ID); the cascode voltages fix the two
    # saturation margins they are derived from.
    RECORNER_RESOLVE = ["M1a_ID", "M5a_ID_over_M1a_ID", "M5a_VDSAT_MARGIN", "M2a_VDSAT_MARGIN"]

    # ---------- (3) UNKNOWNS (node VDS/VSB the DC solver finds) ----------
    UNKNOWNS = [
        Unknown("M1a_VDS", seed=lambda c: c["vdd"] / 3, bound=lambda c: (0.02, c["vdd"])),
        Unknown("M1a_VSB", seed=lambda c: 0, bound=lambda c: (0.02, c["vdd"])),
        Unknown("M2a_VDS", seed=lambda c: c["vdd"] / 4, bound=lambda c: (0.02, c["vdd"])),
        Unknown("M3a_VDS", seed=lambda c: c["vdd"] / 4, bound=lambda c: (0.02, c["vdd"])),
        Unknown("M4a_VDS", seed=lambda c: c["vdd"] / 4, bound=lambda c: (0.02, c["vdd"])),
        Unknown("M5a_VDS", seed=lambda c: c["vdd"] / 4, bound=lambda c: (0.02, c["vdd"])),
        Unknown("M3b_VSB", seed=lambda c: 0, bound=lambda c: (0.02, c["vdd"])),
        Unknown("M4b_VDS", seed=lambda c: c["vdd"] / 4, bound=lambda c: (0.02, c["vdd"])),
        Unknown("M4b_VSB", seed=lambda c: 0, bound=lambda c: (0.02, c["vdd"])),
        Unknown("Mvbn_GMID", seed=lambda c: 12.0, bound=lambda c: (4.0, 30.0)),  # diode replica gm/ID [1]
        Unknown("Mvbp_GMID", seed=lambda c: 12.0, bound=lambda c: (4.0, 30.0)),  # diode replica gm/ID [1]
    ]

    # ---------- (4) SOLVE_POINT ----------
    def solve_point(self, v: State, dev, cond) -> State:
        COUT = cond["cout"]
        VDD = cond["vdd"]
        VIN_CM = cond["vin_cm"]

        M1a = dev.nmos(gmid=v.M1a_GMID, L=v.M1a_L, vds=v.M1a_VDS, vsb=v.M1a_VSB)
        M2a = dev.pmos(gmid=v.M2a_GMID, L=v.M2a_L, vds=v.M2a_VDS, vsb=0.0)
        M5a = dev.nmos(gmid=v.M5a_GMID, L=v.M5a_L, vds=v.M5a_VDS, vsb=0.0)
        M6_VDS = VIN_CM - M1a.vgs
        M6 = dev.nmos(gmid=v.M6_GMID, L=v.M6_L, vds=M6_VDS, vsb=0.0)
        M3a_VSB = M2a.vds_used
        M3a = dev.pmos(gmid=v.M3a_GMID, L=v.M3a_L, vds=v.M3a_VDS, vsb=M3a_VSB)
        M4a_VSB = M5a.vds_used
        M4a = dev.nmos(gmid=v.M4a_GMID, L=v.M4a_L, vds=v.M4a_VDS, vsb=M4a_VSB)
        M3b_VDS = M3a.vds_used
        M3b = dev.pmos(gmid=v.M3a_GMID, L=v.M3a_L, vds=M3b_VDS, vsb=v.M3b_VSB)
        M4b = dev.nmos(gmid=v.M4a_GMID, L=v.M4a_L, vds=v.M4b_VDS, vsb=v.M4b_VSB)
        M5b_VDS = M4a.vgs - M4b.vgs + M5a.vds_used
        M5b = dev.nmos(gmid=v.M5a_GMID, L=v.M5a_L, vds=M5b_VDS, vsb=0.0)
        M2b_VDS = M2a.vds_used + M3a.vgs - M3b.vgs
        M2b = dev.pmos(gmid=v.M2a_GMID, L=v.M2a_L, vds=M2b_VDS, vsb=0.0)
        M1b_VDS = VDD - VIN_CM + M1a.vgs - M2a.vds_used - M3a.vgs + M3b.vgs
        M1b_VSB = VIN_CM - M1a.vgs
        M1b = dev.nmos(gmid=v.M1a_GMID, L=v.M1a_L, vds=M1b_VDS, vsb=M1b_VSB)
        Mvbn = dev.nmos(gmid=v.Mvbn_GMID, L=v.M6_L, vds=M6.vgs, vsb=0.0)
        Mvbp = dev.pmos(gmid=v.Mvbp_GMID, L=v.M2a_L, vds=M2a.vgs, vsb=0.0)

        ID = {}
        ID["M1a"] = v.M1a_ID
        ID["M1b"] = ID["M1a"]
        ID["M2a"] = ID["M1a"] * (v.M5a_ID_over_M1a_ID + 1)
        ID["M3a"] = v.M5a_ID_over_M1a_ID * ID["M1a"]
        ID["M4a"] = v.M5a_ID_over_M1a_ID * ID["M1a"]
        ID["M5a"] = v.M5a_ID_over_M1a_ID * ID["M1a"]
        ID["M6"] = ID["M1a"] + ID["M1b"]
        ID["M2b"] = ID["M2a"]
        ID["M3b"] = ID["M3a"]
        ID["M4b"] = ID["M4a"]
        ID["M5b"] = ID["M5a"]

        W = {}
        W["M1a"] = ID["M1a"] / M1a.jd
        W["M1b"] = W["M1a"]
        W["M2a"] = ID["M2a"] / M2a.jd
        W["M3a"] = ID["M3a"] / M3a.jd
        W["M4a"] = ID["M4a"] / M4a.jd
        W["M5a"] = ID["M5a"] / M5a.jd
        W["M6"] = ID["M6"] / M6.jd
        W["M2b"] = W["M2a"]
        W["M3b"] = W["M3a"]
        W["M4b"] = W["M4a"]
        W["M5b"] = W["M5a"]
        W["Mvbn"] = W["M6"]  # 1:1 bias diode replica
        ID["Mvbn"] = ID["M6"]
        IREF_Mvbn = ID["M6"] * Mvbn.jd / M6.jd
        W["Mvbp"] = W["M2a"]  # 1:1 bias diode replica
        ID["Mvbp"] = ID["M2a"]
        IREF_Mvbp = ID["M2a"] * Mvbp.jd / M2a.jd

        L = {
            "M1a": v.M1a_L,
            "M1b": v.M1a_L,
            "M2a": v.M2a_L,
            "M3a": v.M3a_L,
            "M4a": v.M4a_L,
            "M5a": v.M5a_L,
            "M6": v.M6_L,
            "M2b": v.M2a_L,
            "M3b": v.M3a_L,
            "M4b": v.M4a_L,
            "M5b": v.M5a_L,
            "Mvbn": v.M6_L,
            "Mvbp": v.M2a_L,
        }
        GMID = {
            "M1a": v.M1a_GMID,
            "M1b": v.M1a_GMID,
            "M2a": v.M2a_GMID,
            "M3a": v.M3a_GMID,
            "M4a": v.M4a_GMID,
            "M5a": v.M5a_GMID,
            "M6": v.M6_GMID,
            "M2b": v.M2a_GMID,
            "M3b": v.M3a_GMID,
            "M4b": v.M4a_GMID,
            "M5b": v.M5a_GMID,
            "Mvbn": v.M6_GMID,
            "Mvbp": v.M2a_GMID,
        }

        ptw, ntw = dev.pmos.table_width, dev.nmos.table_width
        ss = {
            "M1a": M1a.small_signal(v.M1a_GMID, ID["M1a"], W["M1a"], ntw, use_gmb=False),
            "M1b": M1b.small_signal(v.M1a_GMID, ID["M1b"], W["M1b"], ntw, use_gmb=False),
            "M2a": M2a.small_signal(v.M2a_GMID, ID["M2a"], W["M2a"], ptw, use_gmb=False),
            "M3a": M3a.small_signal(v.M3a_GMID, ID["M3a"], W["M3a"], ptw, use_gmb=False),
            "M4a": M4a.small_signal(v.M4a_GMID, ID["M4a"], W["M4a"], ntw, use_gmb=False),
            "M5a": M5a.small_signal(v.M5a_GMID, ID["M5a"], W["M5a"], ntw, use_gmb=False),
            "M6": M6.small_signal(v.M6_GMID, ID["M6"], W["M6"], ntw, use_gmb=False),
            "M2b": M2b.small_signal(v.M2a_GMID, ID["M2b"], W["M2b"], ptw, use_gmb=False),
            "M3b": M3b.small_signal(v.M3a_GMID, ID["M3b"], W["M3b"], ptw, use_gmb=False),
            "M4b": M4b.small_signal(v.M4a_GMID, ID["M4b"], W["M4b"], ntw, use_gmb=False),
            "M5b": M5b.small_signal(v.M5a_GMID, ID["M5b"], W["M5b"], ntw, use_gmb=False),
        }

        VBP = VDD - M2a.vgs
        VCASCP = VDD - M2a.vds_used - M3a.vgs
        VCASCN = M4a.vgs + M5a.vds_used
        VBN = M6.vgs
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
            M5a=M5a,
            M6=M6,
            M2b=M2b,
            M3b=M3b,
            M4b=M4b,
            M5b=M5b,
            Mvbn=Mvbn,
            Mvbp=Mvbp,
            IREF_Mvbn=IREF_Mvbn,
            IREF_Mvbp=IREF_Mvbp,
            COUT=COUT,
            VDD=VDD,
            VIN_CM=VIN_CM,
            VBP=VBP,
            VCASCP=VCASCP,
            VCASCN=VCASCN,
            VBN=VBN,
            M5a_ID_over_M1a_ID=v.M5a_ID_over_M1a_ID,
            M5a_VDSAT_MARGIN=v.M5a_VDSAT_MARGIN,
            M2a_VDSAT_MARGIN=v.M2a_VDSAT_MARGIN,
            M1a_VDS=v.M1a_VDS,
            M1a_VSB=v.M1a_VSB,
            M1b_VDS=M1b_VDS,
            M1b_VSB=M1b_VSB,
            M2a_VDS=v.M2a_VDS,
            M3a_VDS=v.M3a_VDS,
            M3a_VSB=M3a_VSB,
            M4a_VDS=v.M4a_VDS,
            M4a_VSB=M4a_VSB,
            M5a_VDS=v.M5a_VDS,
            M6_VDS=M6_VDS,
            M2b_VDS=M2b_VDS,
            M3b_VDS=M3b_VDS,
            M3b_VSB=v.M3b_VSB,
            M4b_VDS=v.M4b_VDS,
            M4b_VSB=v.M4b_VSB,
            M5b_VDS=M5b_VDS,
            Mvbn_GMID=v.Mvbn_GMID,
            Mvbp_GMID=v.Mvbp_GMID,
        )

    # ---------- (5) RESIDUALS (node closures the solver drives to zero) ----------
    def residuals(self, b) -> list:
        return [
            vres(b.M1a_VDS, b.VDD - b.VIN_CM + b.M1a.vgs - b.M2a.vds_used),
            vres(b.M1a_VSB, b.VIN_CM - b.M1a.vgs),
            vres(b.M2a_VDS, b.M2a.vdsat + b.M2a_VDSAT_MARGIN),
            vres(b.M3a_VDS, b.VDD - b.M2a.vds_used - b.M5a.vgs),
            vres(b.M4a_VDS, -b.M5a.vds_used + b.M5a.vgs),
            vres(b.M5a_VDS, b.M5a.vdsat + b.M5a_VDSAT_MARGIN),
            vres(b.M3b_VSB, b.M2a.vds_used + b.M3a.vgs - b.M3b.vgs),
            vres(
                b.M4b_VDS,
                b.VDD
                - b.M2a.vds_used
                - b.M3a.vgs
                - b.M3b.vds_used
                + b.M3b.vgs
                - b.M4a.vgs
                + b.M4b.vgs
                - b.M5a.vds_used,
            ),
            vres(b.M4b_VSB, b.M4a.vgs - b.M4b.vgs + b.M5a.vds_used),
            vres(b.Mvbn.vgs, b.M6.vgs),  # diode replica shares the master's VGS
            vres(b.Mvbp.vgs, b.M2a.vgs),  # diode replica shares the master's VGS
        ]

    # ---------- (7) MULTICORNER CONSERVATION ----------
    def freeze_extra(self, b) -> dict:
        return {"IREF_Mvbn": b.IREF_Mvbn, "IREF_Mvbp": b.IREF_Mvbp, "VCASCN": b.VCASCN, "VCASCP": b.VCASCP}

    def recorner_residuals(self, b, frozen) -> list:
        e = frozen["extra"]
        return [
            rres(b.IREF_Mvbn, e["IREF_Mvbn"], e["IREF_Mvbn"]),
            rres(b.IREF_Mvbp, e["IREF_Mvbp"], e["IREF_Mvbp"]),
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
            + b.L["M5a"] * b.W["M5a"]
            + b.L["M5b"] * b.W["M5b"]
            + b.L["M6"] * b.W["M6"]
            + b.L["Mvbn"] * b.W["Mvbn"]
            + b.L["Mvbp"] * b.W["Mvbp"]
        )
        out["Area"] = Area
        Itotal = b.M5a_ID_over_M1a_ID * b.ID["M1a"] + b.ID["M1a"] + b.ID["M2b"]
        out["Itotal"] = Itotal
        VBP = b.VDD - b.M2a.vgs
        out["VBP"] = VBP
        VCASCP = b.VDD - b.M2a.vds_used - b.M3a.vgs
        out["VCASCP"] = VCASCP
        VCASCN = b.M4a.vgs + b.M5a.vds_used
        out["VCASCN"] = VCASCN
        VBN = b.M6.vgs
        out["VBN"] = VBN
        VOUT_DC = b.VDD - b.M2a.vds_used - b.M3a.vgs - b.M3b.vds_used + b.M3b.vgs
        out["VOUT_DC"] = VOUT_DC
        VOUT_MAX = -b.M3b.vdsat + b.M3b.vgs + b.VCASCP
        out["VOUT_MAX"] = VOUT_MAX
        VOUT_MIN = b.M4b.vdsat - b.M4b.vgs + b.VCASCN
        out["VOUT_MIN"] = VOUT_MIN
        VIN_MAX = min(
            -b.M1a.vdsat + b.M1a.vgs + b.M3a.vgs + b.VCASCP,
            b.M1a.vgs - b.M1b.vdsat + b.M3b.vgs + b.VCASCP,
        )
        out["VIN_MAX"] = VIN_MAX
        VIN_MIN = b.M1a.vgs + b.M6.vdsat
        out["VIN_MIN"] = VIN_MIN
        Output_Swing = VOUT_MAX - VOUT_MIN
        out["Output_Swing"] = Output_Swing
        # saturation margins of the devices whose VDS the fixed cascode biases set: re-solved
        # on every non-reference corner, so they show what the corner does to them
        out["M2a margin"] = b.M2a_VDSAT_MARGIN
        out["M5a margin"] = b.M5a_VDSAT_MARGIN
        out["Margin min"] = min(b.M2a_VDSAT_MARGIN, b.M5a_VDSAT_MARGIN)
        return out

    # ---------- netlist hooks ----------
    def vsource_values(self, ref_op, frozen) -> dict:
        return {"VCASCP": ref_op.VCASCP, "VCASCN": ref_op.VCASCN}

    def mirror_currents(self, ref_op) -> dict:
        return {"VBN": ref_op.IREF_Mvbn, "VBP": ref_op.IREF_Mvbp}

    def netlist_context(self, corner, ref_op=None) -> dict:
        return {"vcm": corner.cond("vin_cm")}


__all__ = ["Circuit", "run", "Knob", "Spec"]
