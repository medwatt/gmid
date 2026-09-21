from __future__ import annotations

import numpy as np

from mosplot.optimizer import (
    CircuitModel,
    Instance,
    Knob,
    Passive,
    State,
    VSource,
    build_ss_model,
    run,
    vres,
)


class Circuit(CircuitModel):
    NAME = "amp"
    PORTS = ["VIN", "VOUT", "vdd", "vss"]
    GROUND = "vss"

    MOSFETS = [
        Instance("M1", "nmos", d="VOUT", g="VIN", s="vss", b="vss"),
    ]

    PASSIVES = [
        Passive("R1", "res", a="vdd", b="VOUT"),
        Passive("COUT", "cap", a="VOUT", b="vss", external=True),
    ]

    VSOURCES = [
        VSource("VDD", p="vdd", n="vss", supply=True),
    ]

    SIGNAL_NODES = {"VIN"}

    KNOBS = [
        Knob("M1_GMID", role="op", sets_width_of="M1"),
        Knob("M1_L", role="geom"),
        Knob("R1", role="geom"),
        Knob("VOUT_Q", role="external"),
    ]

    RECORNER_RESOLVE = ["VOUT_Q"]

    UNKNOWNS = []

    def solve_point(self, v, dev, cond):
        VDD = cond["vdd"]
        VOUT_DC = v.VOUT_Q
        COUT = cond["cout"]

        M1_ID = (VDD - VOUT_DC) / v.R1

        M1 = dev.nmos(gmid=v.M1_GMID, L=v.M1_L, vds=VOUT_DC, vsb=0.0)

        W = {"M1": M1_ID / M1.jd}
        L = {"M1": v.M1_L}
        GMID = {"M1": v.M1_GMID}
        ID = {"M1": M1_ID}

        ntw = dev.nmos.table_width
        ss = {
            "M1": M1.small_signal(v.M1_GMID, ID["M1"], W["M1"], ntw, use_gmb=True),
        }

        VIN_DC = M1.vgs

        return State(
            M1=M1,
            W=W,
            L=L,
            ID=ID,
            GMID=GMID,
            ss=ss,
            VDD=VDD,
            VOUT_DC=VOUT_DC,
            COUT=COUT,
            VIN_DC=VIN_DC,
            R1=v.R1
        )

    def residuals(self, b):
        return []

    def specs(self, b, cond):
        ss_model = build_ss_model(
            self.MOSFETS,
            self.PASSIVES,
            self.VSOURCES,
            b.ss,
            {"R1": b.R1, "COUT": cond["cout"]},
            signal_nodes=self.SIGNAL_NODES,
        )
        ac = ss_model.transfer(inputs={"VIN": 1.0}, output={"VOUT": 1.0})

        Area = b.L["M1"] * b.W["M1"]
        VOUT_MAX = b.VDD
        VOUT_MIN = b.M1.vdsat
        Output_Swing = VOUT_MAX - VOUT_MIN

        return {
            "GBW": ac.ugf(),
            "AC Gain (dB)": 20.0 * np.log10(max(abs(ac.gain()), 1e-300)),
            "PM": ac.phase_margin(),
            "Area": Area,
            "Itotal": b.ID["M1"],
            "VIN_DC": b.VIN_DC,
            "VOUT_DC": b.VOUT_DC,
            "VOUT_MAX": VOUT_MAX,
            "VOUT_MIN": VOUT_MIN,
            "Output_Swing": Output_Swing,
        }

    def freeze_extra(self, b) -> dict:
        return {"VIN_DC": b.VIN_DC}

    def recorner_residuals(self, b, frozen) -> list:
        return [vres(b.VIN_DC, frozen["extra"]["VIN_DC"], 0.05)]

    def passive_values(self, ref_op, frozen=None):
        return {"R1": ref_op.R1}

    def netlist_context(self, corner, ref_op=None) -> dict:
        vcm = ref_op.VIN_DC if ref_op is not None else corner.cond("vdd") / 2
        return {"vcm": vcm, "vout_dc": corner.cond("vout_dc")}


__all__ = ["Circuit", "run", "Knob"]
