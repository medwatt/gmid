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
    PORTS = ["VOUT", "vdd", "vss"]
    GROUND = "vss"

    MOSFETS = [
        Instance("M1", "pmos", d="n01", g="VOUT", s="n02", b="vdd"),
        Instance("M2", "pmos", d="VOUT", g="VOUT", s="vdd", b="vdd"),
        Instance("M3", "nmos", d="n01", g="n01", s="vss", b="vss"),
        Instance("M4", "nmos", d="VOUT", g="n01", s="vss", b="vss"),
    ]
    PASSIVES = [Passive("R1", "res", a="n02", b="vdd")]
    VSOURCES = [VSource("VDD", p="vdd", n="vss", supply=True)]
    SIGNAL_NODES = set()

    KNOBS = [
        Knob("M2_GMID", role="op", sets_width_of="M2"),
        Knob("M3_GMID", role="op", sets_width_of="M3"),
        Knob("L_p", role="geom"),
        Knob("L_n", role="geom"),
        Knob("K", role="geom"),
        Knob("IREF", role="external"),
    ]

    RECORNER_RESOLVE = ["IREF"]

    UNKNOWNS = [
        Unknown("VOUT", seed=lambda c: 0.55 * c["vdd"], bound=lambda c: (0.05, c["vdd"] - 0.05)),
        Unknown("VN", seed=lambda c: 0.40 * c["vdd"], bound=lambda c: (0.05, c["vdd"] - 0.05)),
        Unknown("VR", seed=lambda c: 0.05 * c["vdd"], bound=lambda c: (1e-4, 0.5 * c["vdd"])),
        Unknown("M1_GMID", seed=lambda c: 14.0, bound=(3.0, 24.0)),
        Unknown("M4_GMID", seed=lambda c: 14.0, bound=(3.0, 24.0)),
        Unknown("I_B", seed=lambda c: 20e-6, bound=(1e-6, 1e-3)),
    ]

    def solve_point(self, v: State, dev, cond) -> State:
        VDD = cond["vdd"]
        I_A = v.IREF
        I_B = v.I_B

        M2 = dev.pmos(gmid=v.M2_GMID, L=v.L_p, vds=VDD - v.VOUT, vsb=0.0)
        M3 = dev.nmos(gmid=v.M3_GMID, L=v.L_n, vds=v.VN, vsb=0.0)
        M1 = dev.pmos(gmid=v.M1_GMID, L=v.L_p, vds=VDD - v.VR - v.VN, vsb=v.VR)
        M4 = dev.nmos(gmid=v.M4_GMID, L=v.L_n, vds=v.VOUT, vsb=0.0)

        W = {"M2": I_B / M2.jd, "M3": I_A / M3.jd}
        W["M1"] = v.K * W["M2"]
        W["M4"] = W["M3"]
        R1 = v.VR / I_A

        return State(
            M1=M1,
            M2=M2,
            M3=M3,
            M4=M4,
            W=W,
            L={"M1": v.L_p, "M2": v.L_p, "M3": v.L_n, "M4": v.L_n},
            ID={"M1": I_A, "M2": I_B, "M3": I_A, "M4": I_B},
            GMID={"M1": v.M1_GMID, "M2": v.M2_GMID, "M3": v.M3_GMID, "M4": v.M4_GMID},
            IA=I_A,
            I_B=I_B,
            K=v.K,
            R1=R1,
            VOUT=v.VOUT,
            VN=v.VN,
            VR=v.VR,
            VDD=VDD,
        )

    def residuals(self, b) -> list:
        return [
            vres(b.VDD - b.VOUT, b.M2.vgs),              # M2 diode: VSG2 = VDD - VOUT
            vres(b.VN, b.M3.vgs),                        # M3 diode: VGS3 = VN
            vres(b.VDD - b.VR - b.VOUT, b.M1.vgs),       # M1 VSG set by resistor/gate
            rres(b.IA, b.K * b.W["M2"] * b.M1.jd, b.IA), # M1 carries branch-A current
            vres(b.VN, b.M4.vgs),                        # M4 gate tied to M3 diode gate
            rres(b.I_B, b.W["M4"] * b.M4.jd, b.I_B),     # NMOS mirror output current
        ]

    def freeze_extra(self, b) -> dict:
        return {"R1": b.R1}

    def recorner_residuals(self, b, frozen) -> list:
        return [rres(b.VR, b.IA * frozen["extra"]["R1"], b.VR)]

    def specs(self, b, cond) -> dict:
        gm2 = b.GMID["M2"] * b.I_B
        crit = [b.M1, b.M2, b.M3, b.M4]
        vds_margin = min(pt.vds_used - pt.vdsat for pt in crit)
        return {
            "Area": sum(b.L[d] * b.W[d] for d in b.W),
            "Itotal": b.IA + b.I_B,
            "gm_ref": gm2,
            "gmR": gm2 * b.R1,
            "VDS_margin": vds_margin,
            "I_match": abs(b.I_B - b.IA) / b.IA,
            "VOUT": b.VOUT,
            "VN": b.VN,
            "VR": b.VR,
            "R1": b.R1,
            "IREF": b.IA,
            "I_B": b.I_B,
            "PMOS_GMID": b.GMID["M2"],
            "NMOS_GMID": b.GMID["M3"],
        }

    def passive_values(self, ref_op) -> dict:
        return {"R1": ref_op.R1}


__all__ = ["Circuit", "run", "Knob", "Spec"]
