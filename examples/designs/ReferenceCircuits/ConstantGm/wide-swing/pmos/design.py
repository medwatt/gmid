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
    PORTS = ["vout", "vdd", "vss"]
    GROUND = "vss"

    # ---------- (1) TOPOLOGY ----------
    MOSFETS = [
        Instance("M1", "pmos", d="n4", g="vout", s="n2", b="vdd"),
        Instance("M2", "pmos", d="n5", g="vout", s="vdd", b="vdd"),
        Instance("M3", "pmos", d="n7", g="n6", s="n4", b="vdd"),
        Instance("M4", "pmos", d="vout", g="n6", s="n5", b="vdd"),
        Instance("M5", "nmos", d="n8", g="n7", s="n9", b="vss"),
        Instance("M6", "nmos", d="n6", g="n7", s="n10", b="vss"),
        Instance("M7", "nmos", d="n9", g="n8", s="vss", b="vss"),
        Instance("M8", "nmos", d="n10", g="n8", s="vss", b="vss"),
    ]
    PASSIVES = [
        Passive("R1", "res", a="vdd", b="n2"),
        Passive("Rn", "res", a="n7", b="n8"),
        Passive("Rp", "res", a="vout", b="n6"),
    ]
    VSOURCES = [VSource("VDD", p="vdd", n="vss", supply=True)]
    SIGNAL_NODES = set()

    # ---------- (2) KNOBS ----------
    KNOBS = [
        Knob("M2_GMID", role="op", sets_width_of="M2"),
        Knob("M3_GMID", role="op", sets_width_of="M3"),
        Knob("M5_GMID", role="op", sets_width_of="M5"),
        Knob("M7_GMID", role="op", sets_width_of="M7"),
        Knob("L_p", role="geom"),
        Knob("L_pc", role="geom"),
        Knob("L_nc", role="geom"),
        Knob("L_n", role="geom"),
        Knob("K", role="geom"),
        Knob("IREF", role="external"),
        Knob("Rn_V", role="external"),
        Knob("Rp_V", role="external"),
    ]
    RECORNER_RESOLVE = ["IREF", "Rn_V", "Rp_V"]

    # ---------- (3) UNKNOWNS ----------
    UNKNOWNS = [
        Unknown("V_vout", seed=lambda c: 0.65 * c["vdd"], bound=lambda c: (0.05, c["vdd"])),
        Unknown("V_n2", seed=lambda c: 0.92 * c["vdd"], bound=lambda c: (0.05, c["vdd"])),
        Unknown("V_n4", seed=lambda c: 0.76 * c["vdd"], bound=lambda c: (0.05, c["vdd"])),
        Unknown("V_n5", seed=lambda c: 0.76 * c["vdd"], bound=lambda c: (0.05, c["vdd"])),
        Unknown("V_n6", seed=lambda c: 0.46 * c["vdd"], bound=lambda c: (0.02, c["vdd"])),
        Unknown("V_n7", seed=lambda c: 0.58 * c["vdd"], bound=lambda c: (0.02, c["vdd"])),
        Unknown("V_n8", seed=lambda c: 0.38 * c["vdd"], bound=lambda c: (0.02, c["vdd"])),
        Unknown("V_n9", seed=lambda c: 0.18 * c["vdd"], bound=lambda c: (1e-4, c["vdd"])),
        Unknown("V_n10", seed=lambda c: 0.18 * c["vdd"], bound=lambda c: (1e-4, c["vdd"])),
        Unknown("M1_GMID", seed=lambda c: 14.0, bound=(3.0, 24.0)),
        Unknown("M4_GMID", seed=lambda c: 14.0, bound=(3.0, 24.0)),
        Unknown("M6_GMID", seed=lambda c: 14.0, bound=(3.0, 24.0)),
        Unknown("M8_GMID", seed=lambda c: 14.0, bound=(3.0, 24.0)),
        Unknown("I_R", seed=lambda c: 20e-6, bound=(1e-6, 1e-3)),
    ]

    # ---------- (4) SOLVE_POINT ----------
    def solve_point(self, v: State, dev, cond) -> State:
        VDD = cond["vdd"]
        I_L, I_R, K = v.IREF, v.I_R, v.K

        M2 = dev.pmos(gmid=v.M2_GMID, L=v.L_p, vds=VDD - v.V_n5, vsb=0.0)
        M3 = dev.pmos(gmid=v.M3_GMID, L=v.L_pc, vds=v.V_n4 - v.V_n7, vsb=VDD - v.V_n4)
        M5 = dev.nmos(gmid=v.M5_GMID, L=v.L_nc, vds=v.V_n8 - v.V_n9, vsb=v.V_n9)
        M7 = dev.nmos(gmid=v.M7_GMID, L=v.L_n, vds=v.V_n9, vsb=0.0)

        M1 = dev.pmos(gmid=v.M1_GMID, L=v.L_p, vds=v.V_n2 - v.V_n4, vsb=VDD - v.V_n2)
        M4 = dev.pmos(gmid=v.M4_GMID, L=v.L_pc, vds=v.V_n5 - v.V_vout, vsb=VDD - v.V_n5)
        M6 = dev.nmos(gmid=v.M6_GMID, L=v.L_nc, vds=v.V_n6 - v.V_n10, vsb=v.V_n10)
        M8 = dev.nmos(gmid=v.M8_GMID, L=v.L_n, vds=v.V_n10, vsb=0.0)

        W = {"M2": I_R / M2.jd, "M3": I_L / M3.jd, "M5": I_L / M5.jd, "M7": I_L / M7.jd}
        W["M1"] = K * W["M2"]
        W["M4"] = W["M3"]
        W["M6"] = W["M5"]
        W["M8"] = W["M7"]

        R1 = (VDD - v.V_n2) / I_L
        Rn = (v.V_n7 - v.V_n8) / I_L
        Rp = (v.V_vout - v.V_n6) / I_R

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
                "M1": I_L,
                "M3": I_L,
                "M5": I_L,
                "M7": I_L,
                "M2": I_R,
                "M4": I_R,
                "M6": I_R,
                "M8": I_R,
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
            IL=I_L,
            I_R=I_R,
            K=K,
            R1=R1,
            Rn=Rn,
            Rp=Rp,
            VDD=VDD,
            V_vout=v.V_vout,
            V_n2=v.V_n2,
            V_n4=v.V_n4,
            V_n5=v.V_n5,
            V_n6=v.V_n6,
            V_n7=v.V_n7,
            V_n8=v.V_n8,
            V_n9=v.V_n9,
            V_n10=v.V_n10,
            Rn_V=v.Rn_V,
            Rp_V=v.Rp_V,
        )

    # ---------- (5) RESIDUALS ----------
    def residuals(self, b) -> list:
        return [
            vres(b.M1.vgs, b.V_n2 - b.V_vout),
            vres(b.M2.vgs, b.VDD - b.V_vout),
            vres(b.M3.vgs, b.V_n4 - b.V_n6),
            vres(b.M4.vgs, b.V_n5 - b.V_n6),
            vres(b.M5.vgs, b.V_n7 - b.V_n9),
            vres(b.M6.vgs, b.V_n7 - b.V_n10),
            vres(b.M7.vgs, b.V_n8),
            vres(b.M8.vgs, b.V_n8),
            rres(b.IL, b.W["M1"] * b.M1.jd, b.IL),
            rres(b.I_R, b.W["M4"] * b.M4.jd, b.I_R),
            rres(b.I_R, b.W["M6"] * b.M6.jd, b.I_R),
            rres(b.I_R, b.W["M8"] * b.M8.jd, b.I_R),
            vres(b.V_n7 - b.V_n8, b.Rn_V),
            vres(b.V_vout - b.V_n6, b.Rp_V),
        ]

    # ---------- (6) MULTICORNER CONSERVATION ----------
    def freeze_extra(self, b) -> dict:
        return {"R1": b.R1, "Rn": b.Rn, "Rp": b.Rp}

    def recorner_residuals(self, b, frozen) -> list:
        e = frozen["extra"]
        return [
            rres(b.VDD - b.V_n2, b.IL * e["R1"], b.VDD - b.V_n2),
            rres(b.V_n7 - b.V_n8, b.IL * e["Rn"], b.V_n7 - b.V_n8),
            rres(b.V_vout - b.V_n6, b.I_R * e["Rp"], b.V_vout - b.V_n6),
        ]

    # ---------- (7) SPECS ----------
    def specs(self, b, cond) -> dict:
        gm2 = b.GMID["M2"] * b.I_R
        crit = [b.M1, b.M2, b.M3, b.M4, b.M5, b.M6, b.M7, b.M8]
        vds_margin = min(pt.vds_used - pt.vdsat for pt in crit)
        return {
            "Area": sum(b.L[d] * b.W[d] for d in b.W),
            "Itotal": b.IL + b.I_R,
            "gm_ref": gm2,
            "gmR": gm2 * b.R1,
            "VDS_margin": vds_margin,
            "I_match": abs(b.I_R - b.IL) / b.IL,
            "VOUT": b.V_vout,
            "R1": b.R1,
            "Rn": b.Rn,
            "Rp": b.Rp,
            "IREF": b.IL,
            "I_R": b.I_R,
            "PMOS_GMID": b.GMID["M2"],
            "NMOS_GMID": b.GMID["M7"],
        }

    # ---------- (8) NETLIST HOOKS ----------
    def passive_values(self, ref_op) -> dict:
        return {"R1": ref_op.R1, "Rn": ref_op.Rn, "Rp": ref_op.Rp}


__all__ = ["Circuit", "run", "Knob", "Spec"]
