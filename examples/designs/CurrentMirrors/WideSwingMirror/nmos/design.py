from __future__ import annotations

from mosplot.optimizer import (
    CircuitModel,
    Instance,
    Knob,
    Passive,
    State,
    Unknown,
    VSource,
    build_ss_model,
    run,
    vres,
)


class Circuit(CircuitModel):
    NAME = "wide_swing_mirror"
    PORTS = ["IREF", "VOUT", "vdd", "vss"]
    GROUND = "vss"

    MOSFETS = [
        Instance("M1", "nmos", d="n01", g="IREF", s="gnd", b="gnd"),
        Instance("M2", "nmos", d="n02", g="IREF", s="gnd", b="gnd"),
        Instance("M3", "nmos", d="VOUT", g="vbn", s="n02", b="gnd"),
        Instance("M4", "nmos", d="IREF", g="vbn", s="n01", b="gnd"),
    ]

    PASSIVES = [
        Passive("CL", "cap", a="VOUT", b="gnd", external=True),
    ]

    VSOURCES = [
        VSource("VDD", p="vdd", n="gnd", supply=True),
        VSource("VBN", p="vbn", n="gnd"),
    ]

    SIGNAL_NODES = set()

    KNOBS = [
        Knob("M1_GMID", role="op", sets_width_of="M1"),
        Knob("M1_L", role="geom"),
        Knob("M3_GMID", role="op", sets_width_of="M4"),
        Knob("M3_L", role="geom"),
        Knob("VDSAT_MARGIN", role="external"),   # M1 margin: sets the cascode bias VBN
    ]

    # Multicorner: VBN is emitted as a fixed source, so it is CONSERVED across corners and
    # re-solves the margin it was derived from (the margin then shows what the corner does).
    RECORNER_RESOLVE = ["VDSAT_MARGIN"]

    UNKNOWNS = [
        Unknown("M1_VDS", seed=lambda c: c["vdd"] / 4, bound=lambda c: (0.02, c["vdd"])),
        Unknown("M2_VDS", seed=lambda c: c["vdd"] / 4, bound=lambda c: (0.02, c["vdd"])),
    ]

    def solve_point(self, v, dev, cond):
        VDD = cond["vdd"]
        IREF = cond["iref"]
        K = cond["k"]
        VOUT_DC = cond.get("vout_dc", VDD / 2)
        MARGIN = v.VDSAT_MARGIN

        M1 = dev.nmos(gmid=v.M1_GMID, L=v.M1_L, vds=v.M1_VDS, vsb=0.0)
        M4 = dev.nmos(gmid=v.M3_GMID, L=v.M3_L, vds=M1.vgs - v.M1_VDS, vsb=v.M1_VDS)

        M2 = dev.nmos(gmid=v.M1_GMID, L=v.M1_L, vds=v.M2_VDS, vsb=0.0)
        M3 = dev.nmos(gmid=v.M3_GMID, L=v.M3_L, vds=VOUT_DC - v.M2_VDS, vsb=v.M2_VDS)

        VBN = M4.vgs + v.M1_VDS

        ID = {"M1": IREF, "M2": K * IREF, "M3": K * IREF, "M4": IREF}

        W = {
            "M1": ID["M1"] / M1.jd,
            "M2": ID["M2"] / M2.jd,
            "M3": ID["M3"] / M3.jd,
            "M4": ID["M4"] / M4.jd,
        }

        L = {"M1": v.M1_L, "M2": v.M1_L, "M3": v.M3_L, "M4": v.M3_L}
        GMID = {"M1": v.M1_GMID, "M2": v.M1_GMID, "M3": v.M3_GMID, "M4": v.M3_GMID}

        ntw = dev.nmos.table_width
        ss = {
            "M1": M1.small_signal(v.M1_GMID, ID["M1"], W["M1"], ntw),
            "M2": M2.small_signal(v.M1_GMID, ID["M2"], W["M2"], ntw),
            "M3": M3.small_signal(v.M3_GMID, ID["M3"], W["M3"], ntw),
            "M4": M4.small_signal(v.M3_GMID, ID["M4"], W["M4"], ntw),
        }

        return State(
            M1=M1,
            M2=M2,
            M3=M3,
            M4=M4,
            W=W,
            L=L,
            ID=ID,
            GMID=GMID,
            ss=ss,
            VDD=VDD,
            IREF=IREF,
            K=K,
            VOUT_DC=VOUT_DC,
            VBN=VBN,
            M1_VDS=v.M1_VDS,
            M2_VDS=v.M2_VDS,
            MARGIN=MARGIN,
        )

    def residuals(self, b):
        return [
            vres(b.M1_VDS, b.M1.vdsat + b.MARGIN),
            vres(b.M2_VDS, b.VBN - b.M3.vgs),   # M3's source sits one VGS from the shared bias
        ]

    def specs(self, b, cond):
        ss_model = build_ss_model(
            self.MOSFETS,
            self.PASSIVES,
            self.VSOURCES,
            b.ss,
            {"CL": cond.get("cout", 1e-12)},
        )
        Rout = ss_model.port("VOUT").resistance()

        Area = (
            b.L["M1"] * b.W["M1"]
            + b.L["M2"] * b.W["M2"]
            + b.L["M3"] * b.W["M3"]
            + b.L["M4"] * b.W["M4"]
        )
        Itotal = b.ID["M1"] + b.ID["M2"]
        Vcompliance = b.M2_VDS + b.M3.vdsat

        return {
            "Rout": Rout,
            "Vcompliance": Vcompliance,
            "Iout": b.ID["M2"],
            "Area": Area,
            "Itotal": Itotal,
            # M1 margin (re-solved on other corners: VB is fixed) and M2 margin
            "M1 margin": b.MARGIN,
            "M2 margin": b.M2_VDS - b.M2.vdsat,
            "Margin min": min(b.MARGIN, b.M2_VDS - b.M2.vdsat),
            "VBN": b.VBN,
        }

    def freeze_extra(self, b) -> dict:
        return {"VBN": b.VBN}

    def recorner_residuals(self, b, frozen) -> list:
        return [vres(b.VBN, frozen["extra"]["VBN"], 0.05)]

    def vsource_values(self, ref_op, frozen) -> dict:
        return dict(frozen["extra"])

    def netlist_context(self, corner, ref_op=None) -> dict:
        return {"iref": corner.cond("iref"), "k": corner.cond("k")}


__all__ = ["Circuit", "run", "Knob"]
