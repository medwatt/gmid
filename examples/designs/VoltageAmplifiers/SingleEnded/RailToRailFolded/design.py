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
        Instance("MP1a", "pmos", d="n01", g="VINP", s="ptail", b="vdd"),
        Instance("MP1b", "pmos", d="n03", g="VINN", s="ptail", b="vdd"),
        Instance("MTP", "pmos", d="ptail", g="vbtp", s="vdd", b="vdd"),
        Instance("MN1a", "nmos", d="n05", g="VINP", s="ntail", b="vss"),
        Instance("MN1b", "nmos", d="n06", g="VINN", s="ntail", b="vss"),
        Instance("MTN", "nmos", d="ntail", g="vbtn", s="vss", b="vss"),
        Instance("M2a", "nmos", d="n01", g="vbn", s="vss", b="vss"),
        Instance("M2b", "nmos", d="n03", g="vbn", s="vss", b="vss"),
        Instance("M3a", "nmos", d="n04", g="vcascn", s="n01", b="vss"),
        Instance("M3b", "nmos", d="VOUT", g="vcascn", s="n03", b="vss"),
        Instance("M4a", "pmos", d="n04", g="vcascp", s="n05", b="vdd"),
        Instance("M4b", "pmos", d="VOUT", g="vcascp", s="n06", b="vdd"),
        Instance("M5a", "pmos", d="n05", g="n04", s="vdd", b="vdd"),
        Instance("M5b", "pmos", d="n06", g="n04", s="vdd", b="vdd"),
    ]
    PASSIVES = [
        Passive("COUT", "cap", a="VOUT", b="vss", external=True),
    ]
    VSOURCES = [
        VSource("VDD", p="vdd", n="vss", supply=True),
        VSource("VBTP", p="vbtp", n="vss", mirror="MTP"),
        VSource("VBTN", p="vbtn", n="vss", mirror="MTN"),
        VSource("VBN", p="vbn", n="vss", mirror="M2a"),
        VSource("VCASCN", p="vcascn", n="vss"),
        VSource("VCASCP", p="vcascp", n="vss"),
    ]
    SIGNAL_NODES = {"VINP", "VINN"}

    # ---------- (2) KNOBS ----------
    KNOBS = [
        Knob("MP1a_GMID", role="op", sets_width_of="MP1a"),
        Knob("MN1a_GMID", role="op", sets_width_of="MN1a"),
        Knob("M2a_GMID", role="op", sets_width_of="M2a"),
        Knob("M3a_GMID", role="op", sets_width_of="M3a"),
        Knob("M4a_GMID", role="op", sets_width_of="M4a"),
        Knob("M5a_GMID", role="op", sets_width_of="M5a"),
        Knob("MTP_GMID", role="op", sets_width_of="MTP"),
        Knob("MTN_GMID", role="op", sets_width_of="MTN"),
        Knob("MP1a_L", role="geom"),
        Knob("MN1a_L", role="geom"),
        Knob("M2a_L", role="geom"),
        Knob("M3a_L", role="geom"),
        Knob("M4a_L", role="geom"),
        Knob("M5a_L", role="geom"),
        Knob("MTP_L", role="geom"),
        Knob("MTN_L", role="geom"),
        Knob("M5a_ID_over_MP1a_ID", role="geom"),
        Knob("MN1a_ID_over_MP1a_ID", role="geom"),  # nmos input current / pmos input current
        Knob("MP1a_ID", role="external"),
        Knob("M5a_VDSAT_MARGIN", role="external"),
        Knob("M2a_VDSAT_MARGIN", role="external"),
    ]
    RECORNER_RESOLVE = ["MP1a_ID", "MN1a_ID_over_MP1a_ID", "M5a_ID_over_MP1a_ID",
                        "M2a_VDSAT_MARGIN", "M5a_VDSAT_MARGIN"]


    # ---------- (3) UNKNOWNS ----------
    UNKNOWNS = [
        Unknown("MP1a_VDS", seed=lambda c: c["vdd"] / 3, bound=lambda c: (0.02, c["vdd"])),
        Unknown("M2a_VDS", seed=lambda c: c["vdd"] / 4, bound=lambda c: (0.02, c["vdd"])),
        Unknown("M3a_VDS", seed=lambda c: c["vdd"] / 4, bound=lambda c: (0.02, c["vdd"])),
        Unknown("M4a_VDS", seed=lambda c: c["vdd"] / 4, bound=lambda c: (0.02, c["vdd"])),
        Unknown("M5a_VDS", seed=lambda c: c["vdd"] / 4, bound=lambda c: (0.02, c["vdd"])),
        Unknown("M3b_VSB", seed=lambda c: 0.0, bound=lambda c: (0.02, c["vdd"])),
        Unknown("M4b_VDS", seed=lambda c: c["vdd"] / 4, bound=lambda c: (0.02, c["vdd"])),
        Unknown("M4b_VSB", seed=lambda c: 0.0, bound=lambda c: (0.02, c["vdd"])),
        Unknown("MN1a_VDS", seed=lambda c: c["vdd"] / 3, bound=lambda c: (0.02, c["vdd"])),
        Unknown("MN1_VSB", seed=lambda c: c["vin_cm"] / 3, bound=lambda c: (0.0, c["vdd"])),
        Unknown("Mvbn_GMID", seed=lambda c: 12.0, bound=lambda c: (4.0, 30.0)),
        Unknown("Mvbtp_GMID", seed=lambda c: 12.0, bound=lambda c: (4.0, 30.0)),
        Unknown("Mvbtn_GMID", seed=lambda c: 12.0, bound=lambda c: (4.0, 30.0)),
    ]

    # ---------- (4) SOLVE_POINT ----------
    def solve_point(self, v: State, dev, cond) -> State:
        COUT = cond["cout"]
        VDD = cond["vdd"]
        VIN_CM = cond["vin_cm"]
        VOUT_DC = cond["vout_dc"]

        M2a = dev.nmos(gmid=v.M2a_GMID, L=v.M2a_L, vds=v.M2a_VDS, vsb=0.0)
        M5a = dev.pmos(gmid=v.M5a_GMID, L=v.M5a_L, vds=v.M5a_VDS, vsb=0.0)

        MP1a = dev.pmos(gmid=v.MP1a_GMID, L=v.MP1a_L, vds=v.MP1a_VDS, vsb=0.0)
        MTP_VDS = VDD - VIN_CM - MP1a.vgs
        MTP = dev.pmos(gmid=v.MTP_GMID, L=v.MTP_L, vds=MTP_VDS, vsb=0.0)

        MN1a = dev.nmos(gmid=v.MN1a_GMID, L=v.MN1a_L, vds=v.MN1a_VDS, vsb=v.MN1_VSB)
        MTN_VDS = VIN_CM - MN1a.vgs
        MTN = dev.nmos(gmid=v.MTN_GMID, L=v.MTN_L, vds=MTN_VDS, vsb=0.0)

        M3a_VSB = M2a.vds_used
        M3a = dev.nmos(gmid=v.M3a_GMID, L=v.M3a_L, vds=v.M3a_VDS, vsb=M3a_VSB)
        M4a_VSB = M5a.vds_used
        M4a = dev.pmos(gmid=v.M4a_GMID, L=v.M4a_L, vds=v.M4a_VDS, vsb=M4a_VSB)
        M3b_VDS = M3a.vds_used
        M3b = dev.nmos(gmid=v.M3a_GMID, L=v.M3a_L, vds=M3b_VDS, vsb=v.M3b_VSB)
        M4b = dev.pmos(gmid=v.M4a_GMID, L=v.M4a_L, vds=v.M4b_VDS, vsb=v.M4b_VSB)

        M5b_VDS = M4a.vgs - M4b.vgs + M5a.vds_used
        M5b = dev.pmos(gmid=v.M5a_GMID, L=v.M5a_L, vds=M5b_VDS, vsb=0.0)
        M2b_VDS = M2a.vds_used + M3a.vgs - M3b.vgs
        M2b = dev.nmos(gmid=v.M2a_GMID, L=v.M2a_L, vds=M2b_VDS, vsb=0.0)
        MP1b_VDS = VIN_CM + MP1a.vgs - M2a.vds_used - M3a.vgs + M3b.vgs
        MP1b = dev.pmos(gmid=v.MP1a_GMID, L=v.MP1a_L, vds=MP1b_VDS, vsb=0.0)
        MN1b_VDS = VDD - M5b.vds_used - VIN_CM + MN1a.vgs
        MN1b = dev.nmos(gmid=v.MN1a_GMID, L=v.MN1a_L, vds=MN1b_VDS, vsb=v.MN1_VSB)

        Mvbn = dev.nmos(gmid=v.Mvbn_GMID, L=v.M2a_L, vds=M2a.vgs, vsb=0.0)
        Mvbtp = dev.pmos(gmid=v.Mvbtp_GMID, L=v.MTP_L, vds=MTP.vgs, vsb=0.0)
        Mvbtn = dev.nmos(gmid=v.Mvbtn_GMID, L=v.MTN_L, vds=MTN.vgs, vsb=0.0)

        Ip = v.MP1a_ID
        In = v.MN1a_ID_over_MP1a_ID * Ip
        Ic = v.M5a_ID_over_MP1a_ID * Ip
        I5 = Ic + In

        ID = {}
        ID["MP1a"] = Ip
        ID["MP1b"] = Ip
        ID["MTP"] = 2.0 * Ip
        ID["MN1a"] = In
        ID["MN1b"] = In
        ID["MTN"] = 2.0 * In
        ID["M5a"] = I5
        ID["M5b"] = I5
        ID["M3a"] = Ic
        ID["M3b"] = Ic
        ID["M4a"] = Ic
        ID["M4b"] = Ic
        ID["M2a"] = Ic + Ip
        ID["M2b"] = Ic + Ip

        W = {}
        W["MP1a"] = ID["MP1a"] / MP1a.jd
        W["MP1b"] = W["MP1a"]
        W["MTP"] = ID["MTP"] / MTP.jd
        W["MN1a"] = ID["MN1a"] / MN1a.jd
        W["MN1b"] = W["MN1a"]
        W["MTN"] = ID["MTN"] / MTN.jd
        W["M5a"] = ID["M5a"] / M5a.jd
        W["M5b"] = W["M5a"]
        W["M3a"] = ID["M3a"] / M3a.jd
        W["M3b"] = W["M3a"]
        W["M4a"] = ID["M4a"] / M4a.jd
        W["M4b"] = W["M4a"]
        W["M2a"] = ID["M2a"] / M2a.jd
        W["M2b"] = W["M2a"]

        W["Mvbn"] = W["M2a"]
        ID["Mvbn"] = ID["M2a"]
        IREF_Mvbn = ID["M2a"] * Mvbn.jd / M2a.jd
        W["Mvbtp"] = W["MTP"]
        ID["Mvbtp"] = ID["MTP"]
        IREF_Mvbtp = ID["MTP"] * Mvbtp.jd / MTP.jd
        W["Mvbtn"] = W["MTN"]
        ID["Mvbtn"] = ID["MTN"]
        IREF_Mvbtn = ID["MTN"] * Mvbtn.jd / MTN.jd

        L = {
            "MP1a": v.MP1a_L, "MP1b": v.MP1a_L, "MTP": v.MTP_L,
            "MN1a": v.MN1a_L, "MN1b": v.MN1a_L, "MTN": v.MTN_L,
            "M2a": v.M2a_L, "M2b": v.M2a_L,
            "M3a": v.M3a_L, "M3b": v.M3a_L,
            "M4a": v.M4a_L, "M4b": v.M4a_L,
            "M5a": v.M5a_L, "M5b": v.M5a_L,
            "Mvbn": v.M2a_L, "Mvbtp": v.MTP_L, "Mvbtn": v.MTN_L,
        }
        GMID = {
            "MP1a": v.MP1a_GMID, "MP1b": v.MP1a_GMID, "MTP": v.MTP_GMID,
            "MN1a": v.MN1a_GMID, "MN1b": v.MN1a_GMID, "MTN": v.MTN_GMID,
            "M2a": v.M2a_GMID, "M2b": v.M2a_GMID,
            "M3a": v.M3a_GMID, "M3b": v.M3a_GMID,
            "M4a": v.M4a_GMID, "M4b": v.M4a_GMID,
            "M5a": v.M5a_GMID, "M5b": v.M5a_GMID,
            "Mvbn": v.M2a_GMID, "Mvbtp": v.MTP_GMID, "Mvbtn": v.MTN_GMID,
        }

        ptw, ntw = dev.pmos.table_width, dev.nmos.table_width
        ss = {
            "MP1a": MP1a.small_signal(v.MP1a_GMID, ID["MP1a"], W["MP1a"], ptw, use_gmb=False),
            "MP1b": MP1b.small_signal(v.MP1a_GMID, ID["MP1b"], W["MP1b"], ptw, use_gmb=False),
            "MTP": MTP.small_signal(v.MTP_GMID, ID["MTP"], W["MTP"], ptw, use_gmb=False),
            "MN1a": MN1a.small_signal(v.MN1a_GMID, ID["MN1a"], W["MN1a"], ntw, use_gmb=False),
            "MN1b": MN1b.small_signal(v.MN1a_GMID, ID["MN1b"], W["MN1b"], ntw, use_gmb=False),
            "MTN": MTN.small_signal(v.MTN_GMID, ID["MTN"], W["MTN"], ntw, use_gmb=False),
            "M2a": M2a.small_signal(v.M2a_GMID, ID["M2a"], W["M2a"], ntw, use_gmb=False),
            "M2b": M2b.small_signal(v.M2a_GMID, ID["M2b"], W["M2b"], ntw, use_gmb=False),
            "M3a": M3a.small_signal(v.M3a_GMID, ID["M3a"], W["M3a"], ntw, use_gmb=False),
            "M3b": M3b.small_signal(v.M3a_GMID, ID["M3b"], W["M3b"], ntw, use_gmb=False),
            "M4a": M4a.small_signal(v.M4a_GMID, ID["M4a"], W["M4a"], ptw, use_gmb=False),
            "M4b": M4b.small_signal(v.M4a_GMID, ID["M4b"], W["M4b"], ptw, use_gmb=False),
            "M5a": M5a.small_signal(v.M5a_GMID, ID["M5a"], W["M5a"], ptw, use_gmb=False),
            "M5b": M5b.small_signal(v.M5a_GMID, ID["M5b"], W["M5b"], ptw, use_gmb=False),
        }

        VBN = M2a.vgs
        VBTP = VDD - MTP.vgs
        VBTN = MTN.vgs
        VCASCN = M2a.vds_used + M3a.vgs
        VCASCP = VDD - M4a.vgs - M5a.vds_used

        return State(
            W=W, L=L, ID=ID, GMID=GMID, ss=ss,
            MP1a=MP1a, MP1b=MP1b, MTP=MTP,
            MN1a=MN1a, MN1b=MN1b, MTN=MTN,
            M2a=M2a, M2b=M2b, M3a=M3a, M3b=M3b,
            M4a=M4a, M4b=M4b, M5a=M5a, M5b=M5b,
            Mvbn=Mvbn, Mvbtp=Mvbtp, Mvbtn=Mvbtn,
            IREF_Mvbn=IREF_Mvbn, IREF_Mvbtp=IREF_Mvbtp, IREF_Mvbtn=IREF_Mvbtn,
            COUT=COUT, VDD=VDD, VIN_CM=VIN_CM, VOUT_DC=VOUT_DC,
            VBN=VBN, VBTP=VBTP, VBTN=VBTN, VCASCN=VCASCN, VCASCP=VCASCP,
            M5a_ID_over_MP1a_ID=v.M5a_ID_over_MP1a_ID,
            MN1a_ID_over_MP1a_ID=v.MN1a_ID_over_MP1a_ID,
            M5a_VDSAT_MARGIN=v.M5a_VDSAT_MARGIN,
            M2a_VDSAT_MARGIN=v.M2a_VDSAT_MARGIN,
            MP1a_VDS=v.MP1a_VDS, MP1b_VDS=MP1b_VDS, MTP_VDS=MTP_VDS,
            MN1a_VDS=v.MN1a_VDS, MN1b_VDS=MN1b_VDS, MTN_VDS=MTN_VDS, MN1_VSB=v.MN1_VSB,
            M2a_VDS=v.M2a_VDS, M2b_VDS=M2b_VDS,
            M3a_VDS=v.M3a_VDS, M3a_VSB=M3a_VSB, M3b_VDS=M3b_VDS, M3b_VSB=v.M3b_VSB,
            M4a_VDS=v.M4a_VDS, M4a_VSB=M4a_VSB, M4b_VDS=v.M4b_VDS, M4b_VSB=v.M4b_VSB,
            M5a_VDS=v.M5a_VDS, M5b_VDS=M5b_VDS,
            Mvbn_GMID=v.Mvbn_GMID, Mvbtp_GMID=v.Mvbtp_GMID, Mvbtn_GMID=v.Mvbtn_GMID,
        )

    # ---------- (5) RESIDUALS ----------
    def residuals(self, b) -> list:
        return [
            vres(b.MP1a_VDS, b.VIN_CM + b.MP1a.vgs - b.M2a.vds_used),
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
            vres(b.MN1a_VDS, b.VDD - b.M5a.vds_used - b.VIN_CM + b.MN1a.vgs),
            vres(b.MN1_VSB, b.VIN_CM - b.MN1a.vgs),
            vres(b.Mvbn.vgs, b.M2a.vgs),
            vres(b.Mvbtp.vgs, b.MTP.vgs),
            vres(b.Mvbtn.vgs, b.MTN.vgs),
        ]

    # ---------- (6) MULTICORNER CONSERVATION ----------
    def freeze_extra(self, b) -> dict:
        return {
            "IREF_Mvbn": b.IREF_Mvbn,
            "IREF_Mvbtp": b.IREF_Mvbtp,
            "IREF_Mvbtn": b.IREF_Mvbtn,
            "VCASCN": b.VCASCN,
            "VCASCP": b.VCASCP,
        }

    def recorner_residuals(self, b, frozen) -> list:
        e = frozen["extra"]
        return [
            rres(b.IREF_Mvbn, e["IREF_Mvbn"], e["IREF_Mvbn"]),
            rres(b.IREF_Mvbtp, e["IREF_Mvbtp"], e["IREF_Mvbtp"]),
            rres(b.IREF_Mvbtn, e["IREF_Mvbtn"], e["IREF_Mvbtn"]),
            vres(b.VCASCN, e["VCASCN"], 0.05),
            vres(b.VCASCP, e["VCASCP"], 0.05),
        ]

    # ---------- (7) SPECS ----------
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

        Area = sum(b.L[k] * b.W[k] for k in b.W)
        out["Area"] = Area

        Ip = b.ID["MP1a"]
        In = b.ID["MN1a"]
        Ic = b.ID["M3a"]
        out["Itotal"] = 2.0 * (Ip + In + Ic)

        out["VBN"] = b.VBN
        out["VBTP"] = b.VBTP
        out["VBTN"] = b.VBTN
        out["VCASCN"] = b.VCASCN
        out["VCASCP"] = b.VCASCP

        out["VOUT_DC"] = b.M2a.vds_used + b.M3a.vgs + b.M3b.vds_used - b.M3b.vgs

        VOUT_MAX = -b.M4b.vdsat + b.M4b.vgs + b.VCASCP
        VOUT_MIN = b.M3b.vdsat - b.M3b.vgs + b.VCASCN
        out["VOUT_MAX"] = VOUT_MAX
        out["VOUT_MIN"] = VOUT_MIN
        out["Output_Swing"] = VOUT_MAX - VOUT_MIN

        PMOS_pair_min = max(
            b.MP1a.vdsat - b.MP1a.vgs - b.M3a.vgs + b.VCASCN,
            b.MP1b.vdsat - b.MP1b.vgs - b.M3b.vgs + b.VCASCN,
        )
        NMOS_floor = b.MN1a.vgs + b.MTN.vdsat
        VIN_MIN = max(PMOS_pair_min, NMOS_floor)
        out["VIN_MIN"] = VIN_MIN

        NMOS_pair_max = min(
            b.M4a.vgs - b.MN1a.vdsat + b.MN1a.vgs + b.VCASCP,
            b.M4b.vgs - b.MN1b.vdsat + b.MN1b.vgs + b.VCASCP,
        )
        PMOS_ceiling = b.VDD - b.MP1a.vgs - b.MTP.vdsat
        VIN_MAX = min(PMOS_ceiling, NMOS_pair_max)
        out["VIN_MAX"] = VIN_MAX
        out["ICMR"] = VIN_MAX - VIN_MIN

        out["PMOS_pair_min"] = PMOS_pair_min
        out["NMOS_pair_max"] = NMOS_pair_max
        out["NMOS_floor"] = NMOS_floor
        out["PMOS_ceiling"] = PMOS_ceiling
        out["ICMR_overlap"] = PMOS_ceiling - NMOS_floor
        out["M2a margin"] = b.M2a_VDSAT_MARGIN
        out["M5a margin"] = b.M5a_VDSAT_MARGIN
        out["Margin min"] = min(b.M2a_VDSAT_MARGIN, b.M5a_VDSAT_MARGIN)
        return out

    # ---------- (8) NETLIST HOOKS ----------
    def vsource_values(self, ref_op, frozen) -> dict:
        return {"VCASCN": ref_op.VCASCN, "VCASCP": ref_op.VCASCP}

    def mirror_currents(self, ref_op) -> dict:
        return {
            "VBN": ref_op.IREF_Mvbn,
            "VBTP": ref_op.IREF_Mvbtp,
            "VBTN": ref_op.IREF_Mvbtn,
        }

    def netlist_context(self, corner, ref_op=None) -> dict:
        return {"vcm": corner.cond("vin_cm"), "vout_dc": corner.cond("vout_dc")}


__all__ = ["Circuit", "run", "Knob", "Spec"]
