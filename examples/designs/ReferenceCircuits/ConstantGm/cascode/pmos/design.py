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
    PORTS = ["vbp", "vdd", "vss"]
    GROUND = "vss"

    # ---------- (1) TOPOLOGY ----------
    MOSFETS = [
        Instance("M1", "pmos", d="nA", g="vbp", s="nR", b="vdd"),
        Instance("M2", "pmos", d="vbp", g="vbp", s="vdd", b="vdd"),
        Instance("M3", "pmos", d="ncp", g="ncn", s="nA", b="vdd"),
        Instance("M4", "pmos", d="ncn", g="ncn", s="vbp", b="vdd"),
        Instance("M5", "nmos", d="ncp", g="ncp", s="npm", b="vss"),
        Instance("M6", "nmos", d="ncn", g="ncp", s="nB", b="vss"),
        Instance("M7", "nmos", d="npm", g="npm", s="vss", b="vss"),
        Instance("M8", "nmos", d="nB", g="npm", s="vss", b="vss"),
    ]
    PASSIVES = [Passive("R1", "res", a="vdd", b="nR")]
    VSOURCES = [VSource("VDD", p="vdd", n="vss", supply=True)]
    SIGNAL_NODES = set()

    # ---------- (2) KNOBS ----------
    KNOBS = [
        Knob("M2_GMID", role="op", sets_width_of="M2"),
        Knob("M4_GMID", role="op", sets_width_of="M4"),
        Knob("M5_GMID", role="op", sets_width_of="M5"),
        Knob("M7_GMID", role="op", sets_width_of="M7"),
        Knob("L_p", role="geom"),
        Knob("L_pc", role="geom"),
        Knob("L_nc", role="geom"),
        Knob("L_n", role="geom"),
        Knob("K", role="geom"),
        Knob("IREF", role="external"),
    ]
    RECORNER_RESOLVE = ["IREF"]

    # ---------- (3) UNKNOWNS ----------
    UNKNOWNS = [
        Unknown("V_vbp", seed=lambda c: 0.62 * c["vdd"], bound=lambda c: (0.05, c["vdd"])),
        Unknown("V_ncn", seed=lambda c: 0.38 * c["vdd"], bound=lambda c: (0.05, c["vdd"])),
        Unknown("V_ncp", seed=lambda c: 0.55 * c["vdd"], bound=lambda c: (0.05, c["vdd"])),
        Unknown("V_npm", seed=lambda c: 0.30 * c["vdd"], bound=lambda c: (0.01, c["vdd"])),
        Unknown("V_nA", seed=lambda c: 0.72 * c["vdd"], bound=lambda c: (0.05, c["vdd"])),
        Unknown(
            "V_nR",
            seed=lambda c: 0.95 * c["vdd"],
            bound=lambda c: (0.5 * c["vdd"], c["vdd"] - 1e-4),
        ),
        Unknown("V_nB", seed=lambda c: 0.28 * c["vdd"], bound=lambda c: (0.01, c["vdd"])),
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

        M2 = dev.pmos(gmid=v.M2_GMID, L=v.L_p, vds=VDD - v.V_vbp, vsb=0.0)
        M4 = dev.pmos(gmid=v.M4_GMID, L=v.L_pc, vds=v.V_vbp - v.V_ncn, vsb=VDD - v.V_vbp)
        M5 = dev.nmos(gmid=v.M5_GMID, L=v.L_nc, vds=v.V_ncp - v.V_npm, vsb=v.V_npm)
        M7 = dev.nmos(gmid=v.M7_GMID, L=v.L_n, vds=v.V_npm, vsb=0.0)

        M1 = dev.pmos(gmid=v.M1_GMID, L=v.L_p, vds=v.V_nR - v.V_nA, vsb=VDD - v.V_nR)
        M3 = dev.pmos(gmid=v.M3_GMID, L=v.L_pc, vds=v.V_nA - v.V_ncp, vsb=VDD - v.V_nA)
        M6 = dev.nmos(gmid=v.M6_GMID, L=v.L_nc, vds=v.V_ncn - v.V_nB, vsb=v.V_nB)
        M8 = dev.nmos(gmid=v.M8_GMID, L=v.L_n, vds=v.V_nB, vsb=0.0)

        W = {
            "M2": I_B / M2.jd,
            "M4": I_B / M4.jd,
            "M5": I_A / M5.jd,
            "M7": I_A / M7.jd,
        }
        W["M1"] = K * W["M2"]
        W["M3"] = W["M4"]
        W["M6"] = W["M5"]
        W["M8"] = W["M7"]
        R1 = (VDD - v.V_nR) / I_A

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
                "M1": v.L_p,
                "M2": v.L_p,
                "M3": v.L_pc,
                "M4": v.L_pc,
                "M5": v.L_nc,
                "M6": v.L_nc,
                "M7": v.L_n,
                "M8": v.L_n,
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
            R1=R1,
            VDD=VDD,
            V_vbp=v.V_vbp,
            V_nR=v.V_nR,
            V_nA=v.V_nA,
            V_ncn=v.V_ncn,
            V_ncp=v.V_ncp,
            V_npm=v.V_npm,
            V_nB=v.V_nB,
        )

    # ---------- (5) RESIDUALS ----------
    def residuals(self, b) -> list:
        return [
            vres(b.VDD - b.V_vbp, b.M2.vgs),
            vres(b.V_vbp - b.V_ncn, b.M4.vgs),
            vres(b.V_ncp - b.V_npm, b.M5.vgs),
            vres(b.V_npm, b.M7.vgs),
            vres(b.M1.vgs, b.V_nR - b.V_vbp),
            vres(b.M3.vgs, b.V_nA - b.V_ncn),
            vres(b.M6.vgs, b.V_ncp - b.V_nB),
            vres(b.M8.vgs, b.V_npm),
            rres(b.IA, b.W["M1"] * b.M1.jd, b.IA),
            rres(b.IA, b.W["M3"] * b.M3.jd, b.IA),
            rres(b.I_B, b.W["M6"] * b.M6.jd, b.I_B),
            rres(b.I_B, b.W["M8"] * b.M8.jd, b.I_B),
        ]

    # ---------- (6) MULTICORNER CONSERVATION ----------
    def freeze_extra(self, b) -> dict:
        return {"R1": b.R1}

    def recorner_residuals(self, b, frozen) -> list:
        return [rres(b.VDD - b.V_nR, b.IA * frozen["extra"]["R1"], b.VDD - b.V_nR)]

    # ---------- (7) SPECS ----------
    def specs(self, b, cond) -> dict:
        gm2 = b.GMID["M2"] * b.I_B
        crit = [b.M1, b.M2, b.M3, b.M4, b.M5, b.M6, b.M7, b.M8]
        vds_margin = min(pt.vds_used - pt.vdsat for pt in crit)
        return {
            "Area": sum(b.L[d] * b.W[d] for d in b.W),
            "Itotal": b.IA + b.I_B,
            "gm_ref": gm2,
            "gmR": gm2 * b.R1,
            "VDS_margin": vds_margin,
            "I_match": abs(b.I_B - b.IA) / b.IA,
            "VBP": b.V_vbp,
            "R1": b.R1,
            "IREF": b.IA,
            "I_B": b.I_B,
            "PMOS_GMID": b.GMID["M2"],
            "NMOS_GMID": b.GMID["M7"],
        }

    # ---------- (8) NETLIST HOOKS ----------
    def passive_values(self, ref_op) -> dict:
        return {"R1": ref_op.R1}


__all__ = ["Circuit", "run", "Knob", "Spec"]
