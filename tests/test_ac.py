"""Small-signal AC engine, checked against closed-form circuit theory.

The workhorse fixture is a single gm stage driving R || C:
    DC gain = -gm*R, pole = 1/(2*pi*R*C), UGF ~ gm/(2*pi*C), PM ~ 90 deg.
"""

import numpy as np
import pytest

from mosplot.optimizer import SmallSignalSolver

GM = 1e-3
R = 100e3  # DC gain = 100
C = 1e-12


def make_gm_stage(gm=GM, r=R, c=C) -> SmallSignalSolver:
    ss = SmallSignalSolver()
    ss.add_mos(gm=gm, gds=1.0 / r, d="out", g="in", s="gnd", b="gnd")
    ss.add_cap(c=c, a="out", b="gnd")
    return ss


class TestTransferBasics:
    def test_dc_gain(self):
        ac = make_gm_stage().transfer(inputs={"in": 1.0}, output={"out": 1.0})
        assert ac.gain() == pytest.approx(-GM * R, rel=1e-9)

    def test_gain_scales_with_input_amplitude(self):
        ac = make_gm_stage().transfer(inputs={"in": 0.5}, output={"out": 1.0})
        assert ac.gain() == pytest.approx(-0.5 * GM * R, rel=1e-9)

    def test_response_at_pole_is_3db_down(self):
        ac = make_gm_stage().transfer(inputs={"in": 1.0}, output={"out": 1.0})
        f_pole = 1.0 / (2 * np.pi * R * C)
        h = ac.response(f_pole)
        assert abs(h) == pytest.approx(GM * R / np.sqrt(2.0), rel=1e-6)

    def test_ugf(self):
        ac = make_gm_stage().transfer(inputs={"in": 1.0}, output={"out": 1.0})
        # Exact single-pole UGF: (gm*R/(2*pi*R*C)) * 1/sqrt(1 - 1/(gmR)^2) ~ gm/(2*pi*C)
        assert ac.ugf() == pytest.approx(GM / (2 * np.pi * C), rel=1e-2)

    def test_phase_margin_single_pole(self):
        ac = make_gm_stage().transfer(inputs={"in": 1.0}, output={"out": 1.0})
        pm = ac.phase_margin()
        expected = 180.0 - np.degrees(np.arctan(GM * R * np.sqrt(1 - 1 / (GM * R) ** 2)))
        # Single-pole system: PM slightly above 90 degrees.
        assert pm == pytest.approx(expected, abs=1.0)

    def test_static_system_has_no_ugf(self):
        ss = SmallSignalSolver()
        ss.add_mos(gm=GM, gds=1.0 / R, d="out", g="in", s="gnd", b="gnd")
        ac = ss.transfer(inputs={"in": 1.0}, output={"out": 1.0})
        assert ac.ugf() == float("inf")  # |gain| > 1 at every frequency

    def test_rejection(self):
        ss = make_gm_stage()
        sig = ss.transfer(inputs={"in": 1.0}, output={"out": 1.0})
        weak = ss.transfer(inputs={"in": 0.01}, output={"out": 1.0})
        assert sig.rejection(weak) == pytest.approx(100.0, rel=1e-9)


class TestResistiveDivider:
    def test_two_resistor_divider(self):
        # Resistors are modeled as gds-only MOS stamps (the builder does the same).
        ss = SmallSignalSolver()
        ss.add_mos(gm=0.0, gds=1e-4, d="in", g="gnd", s="mid", b="gnd")  # 10k top
        ss.add_mos(gm=0.0, gds=3e-4, d="mid", g="gnd", s="gnd", b="gnd")  # 3.33k bottom
        ac = ss.transfer(inputs={"in": 1.0}, output={"mid": 1.0})
        assert ac.gain() == pytest.approx(0.25, rel=1e-9)


class TestGmbStamp:
    def test_gmb_adds_to_effective_transconductance(self):
        # Tie bulk to the gate: the stamp behaves like gm + gmb.
        ss = SmallSignalSolver()
        ss.add_mos(gm=GM, gmb=0.2 * GM, gds=1.0 / R, d="out", g="in", s="gnd", b="in")
        ac = ss.transfer(inputs={"in": 1.0}, output={"out": 1.0})
        assert ac.gain() == pytest.approx(-1.2 * GM * R, rel=1e-9)


class TestPortAnalysis:
    # Port analysis has no AC inputs, so every gate must be tied to a driven
    # or grounded node (build_ss_model AC-grounds bias nodes for exactly this
    # reason); a floating gate would make the system singular.
    @staticmethod
    def _grounded_gate_stage() -> SmallSignalSolver:
        ss = SmallSignalSolver()
        ss.add_mos(gm=GM, gds=1.0 / R, d="out", g="gnd", s="gnd", b="gnd")
        ss.add_cap(c=C, a="out", b="gnd")
        return ss

    def test_resistance_of_single_stage_output(self):
        port = self._grounded_gate_stage().port("out")
        assert port.resistance() == pytest.approx(R, rel=1e-9)

    def test_impedance_falls_at_high_frequency(self):
        port = self._grounded_gate_stage().port("out")
        f_pole = 1.0 / (2 * np.pi * R * C)
        z = port.impedance(10 * f_pole)
        assert abs(z) < 0.2 * R

    def test_diode_connected_resistance(self):
        # gm-diode: rout = 1 / (gm + gds).
        ss = SmallSignalSolver()
        ss.add_mos(gm=GM, gds=1.0 / R, d="out", g="out", s="gnd", b="gnd")
        assert ss.port("out").resistance() == pytest.approx(1.0 / (GM + 1.0 / R), rel=1e-9)

    def test_zero_test_current_rejected(self):
        port = make_gm_stage().port("out")
        with pytest.raises(ValueError, match="non-zero"):
            port.impedance(0.0, test_current=0.0)


class TestComposition:
    def _pair(self):
        ss = SmallSignalSolver()
        ss.add_mos(gm=GM, gds=1.0 / R, d="outp", g="in", s="gnd", b="gnd")
        ss.add_mos(gm=GM, gds=1.0 / R, d="outn", g="gnd", s="in", b="gnd")
        a = ss.transfer(inputs={"in": 1.0}, output={"outp": 1.0})
        b = ss.transfer(inputs={"in": 1.0}, output={"outn": 1.0})
        return a, b

    def test_sub_and_add(self):
        a, b = self._pair()
        diff = a - b
        assert diff.gain() == pytest.approx(a.gain() - b.gain(), rel=1e-12)
        tot = a + b
        assert tot.gain() == pytest.approx(a.gain() + b.gain(), rel=1e-12)

    def test_scale_and_negate(self):
        a, _ = self._pair()
        assert (2.0 * a).gain() == pytest.approx(2.0 * a.gain(), rel=1e-12)
        assert (-a).gain() == pytest.approx(-a.gain(), rel=1e-12)

    def test_weighted_output_equals_composition(self):
        ss = SmallSignalSolver()
        ss.add_mos(gm=GM, gds=1.0 / R, d="outp", g="in", s="gnd", b="gnd")
        ss.add_mos(gm=GM, gds=1.0 / R, d="outn", g="gnd", s="in", b="gnd")
        direct = ss.transfer(inputs={"in": 1.0}, output={"outp": 1.0, "outn": -1.0})
        a = ss.transfer(inputs={"in": 1.0}, output={"outp": 1.0})
        b = ss.transfer(inputs={"in": 1.0}, output={"outn": 1.0})
        assert direct.gain() == pytest.approx((a - b).gain(), rel=1e-12)

    def test_incompatible_composition_rejected(self):
        a, _ = self._pair()
        other = make_gm_stage().transfer(inputs={"in": 1.0}, output={"out": 1.0})
        with pytest.raises(ValueError):
            _ = a + other


class TestValidation:
    def test_no_stamps(self):
        ss = SmallSignalSolver()
        with pytest.raises(RuntimeError, match="No stamps"):
            ss.transfer(inputs={"in": 1.0}, output={"out": 1.0})

    def test_unknown_output_node(self):
        ss = make_gm_stage()
        with pytest.raises(ValueError, match="not an unknown node"):
            ss.transfer(inputs={"in": 1.0}, output={"nonexistent": 1.0})

    def test_empty_inputs_rejected(self):
        ss = make_gm_stage()
        with pytest.raises(ValueError, match="At least one"):
            ss.transfer(inputs={}, output={"out": 1.0})

    def test_probe_equals_reference_rejected(self):
        ss = make_gm_stage()
        with pytest.raises(ValueError, match="different"):
            ss.port("out", reference="out")
