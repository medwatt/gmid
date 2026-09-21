"""build_ss_model: topology -> SmallSignalSolver stamping rules."""

import numpy as np
import pytest

from mosplot.optimizer import Instance, Passive, VSource, build_ss_model

GM = 1e-3
R = 100e3


def common_source(**kwargs):
    """One NMOS common-source stage with a resistive load to vdd."""
    mosfets = [Instance("M1", "nmos", d="out", g="vin", s="vss", b="vss")]
    passives = [Passive("RL", "res", a="out", b="vdd")]
    vsources = [VSource("VDD", p="vdd", n="vss", supply=True)]
    ss_params = {"M1": {"gm": GM, "gds": 0.0}}
    passive_params = {"RL": R}
    defaults = dict(signal_nodes={"vin"})
    defaults.update(kwargs)
    return build_ss_model(mosfets, passives, vsources, ss_params, passive_params, **defaults)


class TestAcGrounding:
    def test_supply_rails_are_ac_ground(self):
        ss = common_source()
        ac = ss.transfer(inputs={"vin": 1.0}, output={"out": 1.0})
        assert ac.gain() == pytest.approx(-GM * R, rel=1e-9)

    def test_signal_nodes_stay_live(self):
        # Without signal_nodes the vin node would still be live here (it is
        # not a vsource), but a bias node listed in VSOURCES becomes ground
        # unless excluded by signal_nodes.
        mosfets = [Instance("M1", "nmos", d="out", g="vb", s="vss", b="vss")]
        passives = [Passive("RL", "res", a="out", b="vdd")]
        vsources = [
            VSource("VDD", p="vdd", n="vss", supply=True),
            VSource("VB", p="vb", n="vss"),
        ]
        ss_params = {"M1": {"gm": GM, "gds": 0.0}}

        grounded = build_ss_model(mosfets, passives, vsources, ss_params, {"RL": R})
        # vb was AC-grounded, so driving it as an input has no effect.
        ac0 = grounded.transfer(inputs={"vb": 1.0}, output={"out": 1.0})
        assert ac0.gain() == 0.0

        live = build_ss_model(
            mosfets, passives, vsources, ss_params, {"RL": R}, signal_nodes={"vb"}
        )
        ac = live.transfer(inputs={"vb": 1.0}, output={"out": 1.0})
        assert ac.gain() == pytest.approx(-GM * R, rel=1e-9)


class TestStamping:
    def test_resistor_becomes_gds_stamp(self):
        # Port analysis needs every gate AC-grounded, so bias the gate from a
        # VSource (as a real bias node would be) instead of a signal node.
        mosfets = [Instance("M1", "nmos", d="out", g="vb", s="vss", b="vss")]
        passives = [Passive("RL", "res", a="out", b="vdd")]
        vsources = [
            VSource("VDD", p="vdd", n="vss", supply=True),
            VSource("VB", p="vb", n="vss"),
        ]
        ss = build_ss_model(
            mosfets, passives, vsources, {"M1": {"gm": GM, "gds": 0.0}}, {"RL": R}
        )
        assert ss.port("out").resistance() == pytest.approx(R, rel=1e-9)

    def test_device_caps_loaded(self):
        mosfets = [Instance("M1", "nmos", d="out", g="vin", s="vss", b="vss")]
        passives = [Passive("RL", "res", a="out", b="vdd")]
        vsources = [VSource("VDD", p="vdd", n="vss", supply=True)]
        cdd = 1e-12
        ss_params = {"M1": {"gm": GM, "gds": 0.0, "cdd": cdd}}
        ss = build_ss_model(
            mosfets, passives, vsources, ss_params, {"RL": R}, signal_nodes={"vin"}
        )
        ac = ss.transfer(inputs={"vin": 1.0}, output={"out": 1.0})
        f_pole = 1.0 / (2 * np.pi * R * cdd)
        assert abs(ac.response(f_pole)) == pytest.approx(GM * R / np.sqrt(2), rel=1e-6)

    def test_negative_cap_values_are_rectified(self):
        # Lookup tables may report negative capacitances; the builder abs()es them.
        mosfets = [Instance("M1", "nmos", d="out", g="vin", s="vss", b="vss")]
        vsources = [VSource("VDD", p="vdd", n="vss", supply=True)]
        cdd = 1e-12
        pos = build_ss_model(
            mosfets, [], vsources,
            {"M1": {"gm": GM, "gds": 1.0 / R, "cdd": cdd}},
            {}, signal_nodes={"vin"},
        )
        neg = build_ss_model(
            mosfets, [], vsources,
            {"M1": {"gm": GM, "gds": 1.0 / R, "cdd": -cdd}},
            {}, signal_nodes={"vin"},
        )
        f = 1.0 / (2 * np.pi * R * cdd)
        a = pos.transfer(inputs={"vin": 1.0}, output={"out": 1.0}).response(f)
        b = neg.transfer(inputs={"vin": 1.0}, output={"out": 1.0}).response(f)
        assert abs(a) == pytest.approx(abs(b), rel=1e-12)

    def test_external_passive_is_stamped(self):
        # external=True affects netlist emission only; the AC model includes it.
        mosfets = [Instance("M1", "nmos", d="out", g="vin", s="vss", b="vss")]
        passives = [Passive("CL", "cap", a="out", b="vss", external=True)]
        vsources = [VSource("VDD", p="vdd", n="vss", supply=True)]
        ss = build_ss_model(
            mosfets,
            passives,
            vsources,
            {"M1": {"gm": GM, "gds": 1.0 / R}},
            {"CL": 1e-12},
            signal_nodes={"vin"},
        )
        ac = ss.transfer(inputs={"vin": 1.0}, output={"out": 1.0})
        assert np.isfinite(ac.ugf()) and ac.ugf() > 0

    def test_in_ss_false_device_is_skipped(self):
        # A bias replica with in_ss=False must not load the AC model.
        mosfets = [
            Instance("M1", "nmos", d="out", g="vin", s="vss", b="vss"),
            Instance("Mrep", "nmos", d="out", g="out", s="vss", b="vss", in_ss=False),
        ]
        vsources = [VSource("VDD", p="vdd", n="vss", supply=True)]
        ss = build_ss_model(
            mosfets,
            [],
            vsources,
            {"M1": {"gm": GM, "gds": 1.0 / R}},  # no entry needed for Mrep
            {},
            signal_nodes={"vin"},
        )
        ac = ss.transfer(inputs={"vin": 1.0}, output={"out": 1.0})
        assert ac.gain() == pytest.approx(-GM * R, rel=1e-9)


class TestModelCache:
    def test_same_topology_reuses_compiled_model(self):
        a = common_source()
        b = common_source()
        assert a._model is b._model

    def test_values_change_per_build(self):
        mosfets = [Instance("M1", "nmos", d="out", g="vin", s="vss", b="vss")]
        passives = [Passive("RL", "res", a="out", b="vdd")]
        vsources = [VSource("VDD", p="vdd", n="vss", supply=True)]
        g1 = build_ss_model(
            mosfets, passives, vsources,
            {"M1": {"gm": GM, "gds": 0.0}}, {"RL": R}, signal_nodes={"vin"},
        ).transfer(inputs={"vin": 1.0}, output={"out": 1.0}).gain()
        g2 = build_ss_model(
            mosfets, passives, vsources,
            {"M1": {"gm": 2 * GM, "gds": 0.0}}, {"RL": R}, signal_nodes={"vin"},
        ).transfer(inputs={"vin": 1.0}, output={"out": 1.0}).gain()
        assert g2 == pytest.approx(2 * g1, rel=1e-12)
