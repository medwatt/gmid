from __future__ import annotations

from mosplot.optimizer import (
    CircuitModel,
    Instance,
    Knob,
    Passive,
    Spec,
    State,
    Unknown,
    VSource,
    run,
    rres,
    vres,
)


class Circuit(CircuitModel):
    NAME = "constant_gm"
    PORTS = ["vbn", "vdd", "vss"]
    GROUND = "vss"

    # ---------- (1) TOPOLOGY ----------
    MOSFETS = [
        Instance("M1", "nmos", d="nA", g="vbn", s="nR", b="vss"),
        Instance("M2", "nmos", d="vbn", g="vbn", s="vss", b="vss"),
        Instance("M3", "nmos", d="ncp", g="ncn", s="nA", b="vss"),
        Instance("M4", "nmos", d="ncn", g="ncn", s="vbn", b="vss"),
        Instance("M5", "pmos", d="ncp", g="ncp", s="npm", b="vdd"),
        Instance("M6", "pmos", d="ncn", g="ncp", s="nB", b="vdd"),
        Instance("M7", "pmos", d="npm", g="npm", s="vdd", b="vdd"),
        Instance("M8", "pmos", d="nB", g="npm", s="vdd", b="vdd"),
    ]
    PASSIVES = [Passive("R", "res", a="nR", b="vss")]
    VSOURCES = [VSource("VDD", p="vdd", n="vss", supply=True)]
    SIGNAL_NODES = set()

    # ---------- (2) KNOBS ----------
    KNOBS = [
        Knob("M2_GMID", role="op", sets_width_of="M2"),  # NMOS mirror diode gm/ID
        Knob("M4_GMID", role="op", sets_width_of="M4"),  # NMOS cascode diode gm/ID
        Knob("M5_GMID", role="op", sets_width_of="M5"),  # PMOS cascode diode gm/ID
        Knob("M7_GMID", role="op", sets_width_of="M7"),  # PMOS mirror diode gm/ID
        Knob("L_n", role="geom"),  # M1,M2 length
        Knob("L_nc", role="geom"),  # M3,M4 length
        Knob("L_pc", role="geom"),  # M5,M6 length
        Knob("L_pm", role="geom"),  # M7,M8 length
        Knob("K", role="geom"),  # M1 is K x wider than M2
        Knob("IREF", role="external"),  # left-branch current I_A
    ]
    RECORNER_RESOLVE = ["IREF"]

    # ---------- (3) UNKNOWNS ----------
    UNKNOWNS = [
        Unknown("V_vbn", seed=lambda c: 0.38 * c["vdd"], bound=lambda c: (0.05, c["vdd"])),
        Unknown("V_ncn", seed=lambda c: 0.62 * c["vdd"], bound=lambda c: (0.05, c["vdd"])),
        Unknown("V_ncp", seed=lambda c: 0.45 * c["vdd"], bound=lambda c: (0.05, c["vdd"])),
        Unknown("V_npm", seed=lambda c: 0.70 * c["vdd"], bound=lambda c: (0.05, c["vdd"])),
        Unknown("V_nA", seed=lambda c: 0.28 * c["vdd"], bound=lambda c: (0.01, c["vdd"])),
        Unknown("V_nR", seed=lambda c: 0.05 * c["vdd"], bound=lambda c: (1e-4, 0.5 * c["vdd"])),
        Unknown("V_nB", seed=lambda c: 0.72 * c["vdd"], bound=lambda c: (0.05, c["vdd"])),
        Unknown("M1_GMID", seed=lambda c: 14.0, bound=(3.0, 24.0)),
        Unknown("M3_GMID", seed=lambda c: 14.0, bound=(3.0, 24.0)),
        Unknown("M6_GMID", seed=lambda c: 14.0, bound=(3.0, 24.0)),
        Unknown("M8_GMID", seed=lambda c: 14.0, bound=(3.0, 24.0)),
        Unknown("I_B", seed=lambda c: 20e-6, bound=(1e-6, 1e-3)),
    ]

    # ---------- (4) SOLVE_POINT ----------
    def solve_point(self, v: State, dev, cond) -> State:
        VDD = cond["vdd"]
        I_A, I_B, K = v.IREF, v.I_B, v.K

        M2 = dev.nmos(gmid=v.M2_GMID, L=v.L_n, vds=v.V_vbn, vsb=0.0)
        M4 = dev.nmos(gmid=v.M4_GMID, L=v.L_nc, vds=v.V_ncn - v.V_vbn, vsb=v.V_vbn)
        M5 = dev.pmos(gmid=v.M5_GMID, L=v.L_pc, vds=v.V_npm - v.V_ncp, vsb=VDD - v.V_npm)
        M7 = dev.pmos(gmid=v.M7_GMID, L=v.L_pm, vds=VDD - v.V_npm, vsb=0.0)

        M1 = dev.nmos(gmid=v.M1_GMID, L=v.L_n, vds=v.V_nA - v.V_nR, vsb=v.V_nR)
        M3 = dev.nmos(gmid=v.M3_GMID, L=v.L_nc, vds=v.V_ncp - v.V_nA, vsb=v.V_nA)
        M6 = dev.pmos(gmid=v.M6_GMID, L=v.L_pc, vds=v.V_nB - v.V_ncn, vsb=VDD - v.V_nB)
        M8 = dev.pmos(gmid=v.M8_GMID, L=v.L_pm, vds=VDD - v.V_nB, vsb=0.0)

        W = {
            "M2": I_B / M2.jd,
            "M4": I_B / M4.jd,  # right-branch diodes carry I_B
            "M5": I_A / M5.jd,
            "M7": I_A / M7.jd,  # left-branch  diodes carry I_A
        }
        W["M1"] = K * W["M2"]  # K x M2
        W["M3"] = W["M4"]  # NMOS cascode 1:1
        W["M6"] = W["M5"]  # PMOS cascode 1:1
        W["M8"] = W["M7"]  # PMOS mirror  1:1
        R = v.V_nR / I_A  # resistor drops V_nR at I_A
        return State(
            M1=M1,
            M2=M2,
            M3=M3,
            M4=M4,
            M5=M5,
            M6=M6,
            M7=M7,
            M8=M8,
            W=W,
            L={
                "M1": v.L_n,
                "M2": v.L_n,
                "M3": v.L_nc,
                "M4": v.L_nc,
                "M5": v.L_pc,
                "M6": v.L_pc,
                "M7": v.L_pm,
                "M8": v.L_pm,
            },
            ID={
                "M1": I_A,
                "M3": I_A,
                "M5": I_A,
                "M7": I_A,
                "M2": I_B,
                "M4": I_B,
                "M6": I_B,
                "M8": I_B,
            },
            GMID={
                "M1": v.M1_GMID,
                "M2": v.M2_GMID,
                "M3": v.M3_GMID,
                "M4": v.M4_GMID,
                "M5": v.M5_GMID,
                "M6": v.M6_GMID,
                "M7": v.M7_GMID,
                "M8": v.M8_GMID,
            },
            IA=I_A,
            I_B=I_B,
            K=K,
            R=R,
            VDD=VDD,
            V_vbn=v.V_vbn,
            V_ncn=v.V_ncn,
            V_ncp=v.V_ncp,
            V_npm=v.V_npm,
            V_nA=v.V_nA,
            V_nR=v.V_nR,
            V_nB=v.V_nB,
        )

    # ---------- (5) RESIDUALS ----------
    def residuals(self, b) -> list:
        return [
            vres(b.V_vbn, b.M2.vgs),  # M2
            vres(b.V_ncn - b.V_vbn, b.M4.vgs),  # M4
            vres(b.V_npm - b.V_ncp, b.M5.vgs),  # M5 (VSD = VSG)
            vres(b.VDD - b.V_npm, b.M7.vgs),  # M7
            vres(b.M1.vgs, b.V_vbn - b.V_nR),  # M1 g=vbn s=nR
            vres(b.M3.vgs, b.V_ncn - b.V_nA),  # M3 g=ncn s=nA
            vres(b.M6.vgs, b.V_nB - b.V_ncp),  # M6 VSG = V_nB - V_ncp (g=ncp)
            vres(b.M8.vgs, b.VDD - b.V_npm),  # M8 VSG = VDD - V_npm
            rres(b.IA, b.W["M1"] * b.M1.jd, b.IA),  # M1 carries I_A
            rres(b.IA, b.W["M3"] * b.M3.jd, b.IA),  # M3 carries I_A
            rres(b.I_B, b.W["M6"] * b.M6.jd, b.I_B),  # M6 carries I_B
            rres(b.I_B, b.W["M8"] * b.M8.jd, b.I_B),  # M8 carries I_B (mirror copy)
        ]

    # ---------- (6) MULTICORNER CONSERVATION ----------
    def freeze_extra(self, b) -> dict:
        return {"R": b.R}

    def recorner_residuals(self, b, frozen) -> list:
        return [rres(b.V_nR, b.IA * frozen["extra"]["R"], b.V_nR)]

    # ---------- (7) SPECS ----------
    def specs(self, b, cond) -> dict:
        gm2 = b.GMID["M2"] * b.I_B  # M2 carries I_B
        crit = {"M1": b.M1, "M3": b.M3, "M6": b.M6, "M8": b.M8}
        vds_margin = min(pt.vds_used - pt.vdsat for pt in crit.values())
        return {
            "Area": sum(b.L[d] * b.W[d] for d in b.W),
            "Itotal": b.IA + b.I_B,
            "gm_ref": gm2,
            "gmR": gm2 * b.R,
            "VDS_margin": vds_margin,  # min saturation headroom
            "I_match": abs(b.I_B - b.IA) / b.IA,  # mirror copy error
            "VBN": b.V_vbn,
            "R": b.R,
            "IREF": b.IA,
            "I_B": b.I_B,
            "NMOS_GMID": b.GMID["M2"],
            "PMOS_GMID": b.GMID["M7"],
        }

    # ---------- (8) NETLIST HOOKS ----------
    def passive_values(self, ref_op) -> dict:
        return {"R": ref_op.V_nR / ref_op.IA}


__all__ = ["Circuit", "run", "Knob", "Spec"]
