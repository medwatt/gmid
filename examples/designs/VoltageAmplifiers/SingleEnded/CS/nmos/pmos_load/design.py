from __future__ import annotations

import numpy as np

from mosplot.optimizer import (
    CircuitModel, Instance, Knob, Passive, Spec, State, Unknown, VSource,
    build_ss_model, run, vres, rres,
)

class Circuit(CircuitModel):
    NAME = 'amp'
    PORTS = ['vss', 'vin', 'VOUT', 'vdd']
    GROUND = "vss"

    # ---------- (1) TOPOLOGY ----------
    MOSFETS = [
        Instance('M1', 'nmos', d='VOUT', g='vin', s='vss', b='vss'),
        Instance('M2', 'pmos', d='VOUT', g='vbp', s='vdd', b='vdd'),
    ]
    PASSIVES = [
        Passive('COUT', 'cap', a='VOUT', b='vss', external=True),
    ]
    VSOURCES = [
        VSource('VDD', p='vdd', n='vss', supply=True),
        VSource('VBP', p='vbp', n='vss', mirror='M2'),
        VSource('VIN', p='vin', n='vss', emit=False),
    ]
    SIGNAL_NODES = {'vin'}

    # ---------- (2) KNOBS ----------
    KNOBS = [
        Knob('M1_GMID', role="op", sets_width_of='M1'),
        Knob('M2_GMID', role="op", sets_width_of='M2'),
        Knob('M1_L', role="geom"),
        Knob('M2_L', role="geom"),
        Knob('M1_ID', role="external"),
        Knob('VOUT_Q', role="external"),
    ]
    RECORNER_RESOLVE = ['M1_ID', 'VOUT_Q']

    # ---------- (3) UNKNOWNS ----------
    UNKNOWNS = [
        Unknown('Mvbp_GMID', seed=lambda c: 12.0, bound=lambda c: (4.0, 30.0)),
    ]

    # ---------- (4) SOLVE_POINT ----------
    def solve_point(self, v: State, dev, cond) -> State:
        COUT = cond['cout']
        VDD = cond['vdd']
        VOUT_DC = v.VOUT_Q

        M1_VDS = VOUT_DC
        M2_VDS = VDD - VOUT_DC

        M1 = dev.nmos(gmid=v.M1_GMID, L=v.M1_L, vds=M1_VDS, vsb=0.0)
        M2 = dev.pmos(gmid=v.M2_GMID, L=v.M2_L, vds=M2_VDS, vsb=0.0)
        Mvbp = dev.pmos(gmid=v.Mvbp_GMID, L=v.M2_L, vds=M2.vgs, vsb=0.0)

        ID = {}
        ID['M1'] = v.M1_ID
        ID['M2'] = ID['M1']

        W = {}
        W['M1'] = ID['M1'] / M1.jd
        W['M2'] = ID['M2'] / M2.jd
        W['Mvbp'] = W['M2']            # 1:1 bias diode replica
        ID['Mvbp'] = ID['M2']
        IREF_Mvbp = ID['M2'] * Mvbp.jd / M2.jd

        L = {
            'M1': v.M1_L,
            'M2': v.M2_L,
            'Mvbp': v.M2_L,
        }
        GMID = {
            'M1': v.M1_GMID,
            'M2': v.M2_GMID,
            'Mvbp': v.Mvbp_GMID,
        }

        ptw, ntw = dev.pmos.table_width, dev.nmos.table_width
        ss = {
            'M1': M1.small_signal(v.M1_GMID, ID['M1'], W['M1'], ntw, use_gmb=False),
            'M2': M2.small_signal(v.M2_GMID, ID['M2'], W['M2'], ptw, use_gmb=False),
        }

        VBP = VDD - M2.vgs
        return State(W=W, L=L, ID=ID, GMID=GMID, ss=ss, M1=M1, M2=M2, Mvbp=Mvbp, IREF_Mvbp=IREF_Mvbp, COUT=COUT, VDD=VDD, VOUT_DC=VOUT_DC, VBP=VBP, M1_VDS=M1_VDS, M2_VDS=M2_VDS, Mvbp_GMID=v.Mvbp_GMID)

    # ---------- (5) RESIDUALS ----------
    def residuals(self, b) -> list:
        return [
            vres(b.Mvbp.vgs, b.M2.vgs),   # diode replica shares the master's VGS
        ]

    # ---------- (6) MULTICORNER CONSERVATION ----------
    def freeze_extra(self, b) -> dict:
        return {'IREF_Mvbp': b.IREF_Mvbp, 'VIN': b.M1.vgs}

    def recorner_residuals(self, b, frozen) -> list:
        e = frozen['extra']
        return [rres(b.IREF_Mvbp, e['IREF_Mvbp'], e['IREF_Mvbp']),
                vres(b.M1.vgs, e['VIN'], 0.05)]

    # ---------- (7) SPECS ----------
    def specs(self, b, cond) -> dict:
        ss = build_ss_model(self.MOSFETS, self.PASSIVES, self.VSOURCES, b.ss,
                            {'COUT': cond['cout']}, signal_nodes=self.SIGNAL_NODES)
        ac = ss.transfer(inputs={'vin': 1.0}, output={'VOUT': 1.0})
        out = {
            "GBW": ac.ugf(),
            "AC Gain (dB)": 20.0 * np.log10(max(abs(ac.gain()), 1e-300)),
            "PM": ac.phase_margin(),
        }
        out['VOUT_DC'] = b.VOUT_DC + cond['vin_cm'] - b.M1.vgs
        Area = b.L['M1']*b.W['M1'] + b.L['M2']*b.W['M2'] + b.L['Mvbp']*b.W['Mvbp']
        out['Area'] = Area
        Itotal = b.ID['M1']
        out['Itotal'] = Itotal
        VBP = b.VDD - b.M2.vgs
        out['VBP'] = VBP
        VIN_DC = b.M1.vgs
        out['VIN_DC'] = VIN_DC
        VOUT_MAX = -b.M2.vdsat + b.VDD
        out['VOUT_MAX'] = VOUT_MAX
        VOUT_MIN = b.M1.vdsat
        out['VOUT_MIN'] = VOUT_MIN
        Output_Swing = VOUT_MAX - VOUT_MIN
        out['Output_Swing'] = Output_Swing
        return out

    # ---------- (8) NETLIST HOOKS ----------
    def mirror_currents(self, ref_op) -> dict:
        return {'VBP': ref_op.IREF_Mvbp}

    def netlist_context(self, corner, ref_op=None) -> dict:
        return {"vcm": corner.cond("vin_cm"), "vout_dc": corner.cond("vout_dc")}

__all__ = ["Circuit", "run", "Knob", "Spec"]
