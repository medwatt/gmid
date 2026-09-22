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
    NAME = "constant_gm_wideswing"
    PORTS = ["vout", "vdd", "vss"]
    GROUND = "vss"

    # ---------- (1) TOPOLOGY ----------
    MOSFETS = [
        Instance("M1", "nmos", d="nLn", g="vout", s="nr", b="vss"),
        Instance("M2", "nmos", d="nRn", g="vout", s="vss", b="vss"),
        Instance("M3", "nmos", d="ncp", g="ncn", s="nLn", b="vss"),
        Instance("M4", "nmos", d="vout", g="ncn", s="nRn", b="vss"),
        Instance("M5", "pmos", d="npm", g="ncp", s="nLp", b="vdd"),
        Instance("M6", "pmos", d="ncn", g="ncp", s="nRp", b="vdd"),
        Instance("M7", "pmos", d="nLp", g="npm", s="vdd", b="vdd"),
        Instance("M8", "pmos", d="nRp", g="npm", s="vdd", b="vdd"),
    ]
    PASSIVES = [
        Passive("R1", "res", a="nr", b="vss"),
        Passive("Rp", "res", a="npm", b="ncp"),
        Passive("Rn", "res", a="ncn", b="vout"),
    ]
    VSOURCES = [VSource("VDD", p="vdd", n="vss", supply=True)]
    SIGNAL_NODES = set()

    # ---------- (2) KNOBS ----------
    KNOBS = [
        Knob("M2_GMID", role="op", sets_width_of="M2"),  # NMOS core (1x) reference
        Knob("M3_GMID", role="op", sets_width_of="M3"),  # NMOS cascode reference
        Knob("M5_GMID", role="op", sets_width_of="M5"),  # PMOS cascode reference
        Knob("M7_GMID", role="op", sets_width_of="M7"),  # PMOS mirror reference
        Knob("L_n", role="geom"),
        Knob("L_nc", role="geom"),
        Knob("L_pc", role="geom"),
        Knob("L_pm", role="geom"),
        Knob("K", role="geom"),  # M1 is K x wider than M2
        Knob("IREF", role="external"),  # left-branch current I_L
        Knob("Rp_V", role="external"),  # wide-swing drop across Rp (I_L*Rp)
        Knob("Rn_V", role="external"),  # wide-swing drop across Rn (I_R*Rn)
    ]
    RECORNER_RESOLVE = ["IREF", "Rp_V", "Rn_V"]

    # ---------- (3) UNKNOWNS ----------
    UNKNOWNS = [
        Unknown("V_vout", seed=lambda c: 0.28 * c["vdd"], bound=lambda c: (0.05, c["vdd"])),
        Unknown("V_nr", seed=lambda c: 0.05 * c["vdd"], bound=lambda c: (1e-4, 0.6 * c["vdd"])),
        Unknown("V_nLn", seed=lambda c: 0.18 * c["vdd"], bound=lambda c: (0.02, c["vdd"])),
        Unknown("V_nRn", seed=lambda c: 0.18 * c["vdd"], bound=lambda c: (0.02, c["vdd"])),
        Unknown("V_ncn", seed=lambda c: 0.45 * c["vdd"], bound=lambda c: (0.05, c["vdd"])),
        Unknown("V_ncp", seed=lambda c: 0.60 * c["vdd"], bound=lambda c: (0.05, c["vdd"])),
        Unknown("V_npm", seed=lambda c: 0.72 * c["vdd"], bound=lambda c: (0.05, c["vdd"])),
        Unknown("V_nLp", seed=lambda c: 0.82 * c["vdd"], bound=lambda c: (0.05, c["vdd"])),
        Unknown("V_nRp", seed=lambda c: 0.82 * c["vdd"], bound=lambda c: (0.05, c["vdd"])),
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

        M2 = dev.nmos(gmid=v.M2_GMID, L=v.L_n, vds=v.V_nRn, vsb=0.0)
        M3 = dev.nmos(gmid=v.M3_GMID, L=v.L_nc, vds=v.V_ncp - v.V_nLn, vsb=v.V_nLn)
        M5 = dev.pmos(gmid=v.M5_GMID, L=v.L_pc, vds=v.V_nLp - v.V_npm, vsb=VDD - v.V_nLp)
        M7 = dev.pmos(gmid=v.M7_GMID, L=v.L_pm, vds=VDD - v.V_nLp, vsb=0.0)

        M1 = dev.nmos(gmid=v.M1_GMID, L=v.L_n, vds=v.V_nLn - v.V_nr, vsb=v.V_nr)
        M4 = dev.nmos(gmid=v.M4_GMID, L=v.L_nc, vds=v.V_vout - v.V_nRn, vsb=v.V_nRn)
        M6 = dev.pmos(gmid=v.M6_GMID, L=v.L_pc, vds=v.V_nRp - v.V_ncn, vsb=VDD - v.V_nRp)
        M8 = dev.pmos(gmid=v.M8_GMID, L=v.L_pm, vds=VDD - v.V_nRp, vsb=0.0)

        W = {"M2": I_R / M2.jd, "M3": I_L / M3.jd, "M5": I_L / M5.jd, "M7": I_L / M7.jd}
        W["M1"] = K * W["M2"]  # K x M2
        W["M4"] = W["M3"]  # NMOS cascode 1:1
        W["M6"] = W["M5"]  # PMOS cascode 1:1
        W["M8"] = W["M7"]  # PMOS mirror  1:1
        R1 = v.V_nr / I_L
        Rp = (v.V_npm - v.V_ncp) / I_L
        Rn = (v.V_ncn - v.V_vout) / I_R
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
            Rp=Rp,
            Rn=Rn,
            VDD=VDD,
            V_vout=v.V_vout,
            V_nr=v.V_nr,
            V_nLn=v.V_nLn,
            V_nRn=v.V_nRn,
            V_ncn=v.V_ncn,
            V_ncp=v.V_ncp,
            V_npm=v.V_npm,
            V_nLp=v.V_nLp,
            V_nRp=v.V_nRp,
            Rp_V=v.Rp_V,
            Rn_V=v.Rn_V,
        )

    # ---------- (5) RESIDUALS ----------
    def residuals(self, b) -> list:
        return [
            vres(b.M1.vgs, b.V_vout - b.V_nr),  # M1 g=vout s=nr
            vres(b.M2.vgs, b.V_vout),  # M2 g=vout s=gnd
            vres(b.M3.vgs, b.V_ncn - b.V_nLn),  # M3 g=ncn s=nLn
            vres(b.M4.vgs, b.V_ncn - b.V_nRn),  # M4 g=ncn s=nRn
            vres(b.M5.vgs, b.V_nLp - b.V_ncp),  # M5 VSG = s(nLp) - g(ncp)
            vres(b.M6.vgs, b.V_nRp - b.V_ncp),  # M6 VSG = s(nRp) - g(ncp)
            vres(b.M7.vgs, b.VDD - b.V_npm),  # M7 VSG = VDD - g(npm)
            vres(b.M8.vgs, b.VDD - b.V_npm),  # M8 VSG = VDD - g(npm)
            rres(b.IL, b.W["M1"] * b.M1.jd, b.IL),  # M1 carries I_L
            rres(b.I_R, b.W["M4"] * b.M4.jd, b.I_R),  # M4 carries I_R
            rres(b.I_R, b.W["M6"] * b.M6.jd, b.I_R),  # M6 carries I_R
            rres(b.I_R, b.W["M8"] * b.M8.jd, b.I_R),  # M8 carries I_R (PMOS mirror copy)
            vres(b.V_npm - b.V_ncp, b.Rp_V),  # Rp drop (left, I_L)
            vres(b.V_ncn - b.V_vout, b.Rn_V),  # Rn drop (right, I_R)
        ]

    # ---------- (6) MULTICORNER CONSERVATION ----------
    def freeze_extra(self, b) -> dict:
        return {"R1": b.R1, "Rp": b.Rp, "Rn": b.Rn}

    def recorner_residuals(self, b, frozen) -> list:
        e = frozen["extra"]
        return [
            rres(b.V_nr, b.IL * e["R1"], b.V_nr),
            rres(b.V_npm - b.V_ncp, b.IL * e["Rp"], b.V_npm - b.V_ncp),
            rres(b.V_ncn - b.V_vout, b.I_R * e["Rn"], b.V_ncn - b.V_vout),
        ]

    # ---------- (7) SPECS ----------
    def specs(self, b, cond) -> dict:
        gm2 = b.GMID["M2"] * b.I_R  # M2 carries I_R
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
            "Rp": b.Rp,
            "Rn": b.Rn,
            "IREF": b.IL,
            "I_R": b.I_R,
            "NMOS_GMID": b.GMID["M2"],
            "PMOS_GMID": b.GMID["M7"],
        }

    # ---------- (8) NETLIST HOOKS ----------
    def passive_values(self, ref_op) -> dict:
        return {"R1": ref_op.R1, "Rp": ref_op.Rp, "Rn": ref_op.Rn}


__all__ = ["Circuit", "run", "Knob", "Spec"]
