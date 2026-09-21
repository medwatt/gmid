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
    PORTS = ["VBN", "VBP2", "vdd", "vss"]
    GROUND = "vss"

    # ============================================================= (1) TOPOLOGY
    MOSFETS = [
        Instance("M1", "nmos", d="VBP2", g="VBN", s="nR", b="vss"),
        Instance("M2", "nmos", d="VBN", g="VBN", s="vss", b="vss"),
        Instance("M3", "pmos", d="VBP2", g="VBP2", s="vdd", b="vdd"),
        Instance("M4", "pmos", d="VBN", g="VBP2", s="vdd", b="vdd"),
    ]
    PASSIVES = [Passive("R", "res", a="nR", b="vss")]
    VSOURCES = [VSource("VDD", p="vdd", n="vss", supply=True)]
    SIGNAL_NODES = set()

    # ============================================================= (2) KNOBS
    KNOBS = [
        Knob("M2_GMID", role="op", sets_width_of="M2"),  # NMOS gm/ID (sets W2)
        Knob("M3_GMID", role="op", sets_width_of="M3"),  # PMOS gm/ID (sets W3)
        Knob("L_n", role="geom"),  # NMOS length (M1 and M2)
        Knob("L_p", role="geom"),  # PMOS length (M3 and M4)
        Knob("K", role="geom"),    # M1 is K x wider than M2
        Knob("IREF", role="external"),  # self-biased current I_A (left branch)
    ]
    # The conserved quantity across corners is the physical resistor R (a DERIVED value),
    # NOT a current source -- so IREF re-solves per corner.
    RECORNER_RESOLVE = ["IREF"]

    # ============================================================= (3) UNKNOWNS
    UNKNOWNS = [
        Unknown("VBN", seed=lambda c: 0.45 * c["vdd"], bound=lambda c: (0.05, c["vdd"])),
        Unknown("VBP2", seed=lambda c: 0.60 * c["vdd"], bound=lambda c: (0.05, c["vdd"])),
        Unknown("VnR", seed=lambda c: 0.05, bound=lambda c: (1e-4, 0.5 * c["vdd"])),
        Unknown("M1_GMID", seed=lambda c: 14.0, bound=(3.0, 24.0)),
        Unknown("M4_GMID", seed=lambda c: 14.0, bound=(3.0, 24.0)),
        Unknown("I_B", seed=lambda c: 20e-6, bound=(1e-6, 1e-3)),
    ]

    # ============================================================= (4) SOLVE_POINT
    def solve_point(self, v: State, dev, cond) -> State:
        VDD = cond["vdd"]
        K = v.K
        I_A = v.IREF # branch A: M3 (diode) -> M1 (K x, degenerated) -> R -> vss
        I_B = v.I_B  # branch B: M4 (mirror) -> M2 (diode) -> vss
        M2 = dev.nmos(gmid=v.M2_GMID, L=v.L_n, vds=v.VBN, vsb=0.0)            # diode  (carries I_B)
        M3 = dev.pmos(gmid=v.M3_GMID, L=v.L_p, vds=VDD - v.VBP2, vsb=0.0)     # diode  (carries I_A)
        M1 = dev.nmos(gmid=v.M1_GMID, L=v.L_n, vds=v.VBP2 - v.VnR, vsb=v.VnR) # body effect (I_A)
        M4 = dev.pmos(gmid=v.M4_GMID, L=v.L_p, vds=VDD - v.VBN, vsb=0.0)      # mirror of M3 (I_B)

        # widths from the reference devices of each branch (diode currents set the widths):
        W = {"M2": I_B / M2.jd, "M3": I_A / M3.jd}
        W["M1"] = K * W["M2"]  # geometry: K x wider than M2
        W["M4"] = W["M3"]  # geometry: 1:1 PMOS mirror
        R = v.VnR / I_A  # resistor drops VnR at the branch-A current
        return State(
            M1=M1,
            M2=M2,
            M3=M3,
            M4=M4,
            W=W,
            L={"M1": v.L_n, "M2": v.L_n, "M3": v.L_p, "M4": v.L_p},
            ID={"M1": I_A, "M2": I_B, "M3": I_A, "M4": I_B},
            GMID={"M1": v.M1_GMID, "M2": v.M2_GMID, "M3": v.M3_GMID, "M4": v.M4_GMID},
            IA=I_A,
            I_B=I_B,
            K=K,
            R=R,
            VBN=v.VBN,
            VBP2=v.VBP2,
            VnR=v.VnR,
            VDD=VDD,
        )

    # ============================================================= (5) RESIDUALS
    def residuals(self, b) -> list:
        return [
            vres(b.VBN, b.M2.vgs),                       # M2 diode: VBN = VGS2
            vres(b.VDD - b.VBP2, b.M3.vgs),              # M3 diode: VDD-VBP2 = VSG3
            vres(b.VBN - b.VnR, b.M1.vgs),               # M1 VGS set by nodes
            rres(b.IA, b.K * b.W["M2"] * b.M1.jd, b.IA), # NMOS copy: M1 (K x W2) carries I_A
            vres(b.VDD - b.VBP2, b.M4.vgs),              # M4 VSG = M3 VSG (gates tied at VBP2)
            rres(
                b.I_B, b.W["M4"] * b.M4.jd, b.I_B
            ),  # PMOS copy: M4 (W4=W3) delivers I_B #   at its own VDS -> finite-lambda
        ]

    # ====================================================== (7) MULTICORNER CONSERVATION
    # The physical resistor R is what is held across corners (a derived value). IREF re-solves
    # so that VnR = IREF * R still holds on each corner's LUT -- the self-bias loop closes.
    def freeze_extra(self, b) -> dict:
        return {"R": b.R}

    def recorner_residuals(self, b, frozen) -> list:
        return [rres(b.VnR, b.IA * frozen["extra"]["R"], b.VnR)]

    # ============================================================= (6) SPECS
    def specs(self, b, cond) -> dict:
        gm2 = b.GMID["M2"] * b.I_B  # M2 carries the branch-B current
        # keep every device saturated (M2/M3 are diodes; M1/M4 are the ones that can go triode):
        crit = [b.M1, b.M2, b.M3, b.M4]
        vds_margin = min(pt.vds_used - pt.vdsat for pt in crit)
        return {
            "Area": sum(b.L[d] * b.W[d] for d in b.W),
            "Itotal": b.IA + b.I_B,
            "gm_ref": gm2,                       # reference transconductance (gm of M2)
            "gmR": gm2 * b.R,                    # constant-gm figure of merit (gm2 * R)
            "VDS_margin": vds_margin,            # min saturation headroom (-> keep >0)
            "I_match": abs(b.I_B - b.IA) / b.IA, # PMOS mirror copy error (supply sensitivity)
            "VBN": b.VBN,
            "VBP2": b.VBP2,
            "R": b.R,
            "IREF": b.IA,
            "I_B": b.I_B,
            "NMOS_GMID": b.GMID["M2"],
            "PMOS_GMID": b.GMID["M3"],
        }

    # ---- netlist hook: R = VnR / IREF (both are solved/knob values on the reference op) ----
    def passive_values(self, ref_op) -> dict:
        return {"R": ref_op.VnR / ref_op.IA}


__all__ = ["Circuit", "run", "Knob", "Spec"]
