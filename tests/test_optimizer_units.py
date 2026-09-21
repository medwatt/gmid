"""Pure-logic units of the optimizer: cost, knob handling, worst-case merge."""

import numpy as np
import pytest

from mosplot.optimizer import CircuitModel, Knob, Spec
from mosplot.optimizer.evaluate import CornerResult, worst_case
from mosplot.optimizer.optimizer import Optimizer, _resolve_knobs, _softplus


def make_optimizer(parameters, target_specs=None):
    class Toy(CircuitModel):
        KNOBS = [Knob(p.name, role=p.role, sets_width_of=p.sets_width_of) for p in parameters]

    return Optimizer(Toy(), parameters, target_specs or {}, corners=[])


class TestSoftplus:
    def test_smooth_near_zero(self):
        assert _softplus(0.0) == pytest.approx(np.log(2) / 20.0)

    def test_linear_for_large_x(self):
        assert _softplus(5.0) == 5.0

    def test_small_for_negative_x(self):
        assert _softplus(-1.0) < 1e-8


class TestComputeCost:
    def _opt(self, specs):
        return make_optimizer([Knob("x", (0.0, 1.0))], specs)

    def test_max_spec_satisfied_is_cheap(self):
        opt = self._opt({"GBW": Spec(10e6, "max", 1.0)})
        good = opt.compute_cost({"GBW": 20e6})
        bad = opt.compute_cost({"GBW": 1e6})
        assert good < 1e-6
        assert bad > good

    def test_max_spec_cost_grows_with_violation(self):
        opt = self._opt({"GBW": Spec(10e6, "max", 1.0)})
        assert opt.compute_cost({"GBW": 1e6}) > opt.compute_cost({"GBW": 5e6})

    def test_min_spec_keeps_pressure_when_satisfied(self):
        opt = self._opt({"I": Spec(10e-6, "min", 1.0)})
        at_half = opt.compute_cost({"I": 5e-6})
        at_quarter = opt.compute_cost({"I": 2.5e-6})
        assert 0 < at_quarter < at_half  # secondary term keeps driving down

    def test_eq_spec_two_sided(self):
        opt = self._opt({"V": Spec(0.6, "eq", 1.0)})
        assert opt.compute_cost({"V": 0.6}) == pytest.approx(0.0)
        assert opt.compute_cost({"V": 0.5}) > 0
        assert opt.compute_cost({"V": 0.7}) > 0

    def test_missing_spec_is_heavily_penalised(self):
        opt = self._opt({"V": Spec(0.6, "eq", 1.0)})
        assert opt.compute_cost({}) >= 1e6

    def test_weight_scales_cost(self):
        light = self._opt({"GBW": Spec(10e6, "max", 1.0)}).compute_cost({"GBW": 1e6})
        heavy = self._opt({"GBW": Spec(10e6, "max", 4.0)}).compute_cost({"GBW": 1e6})
        assert heavy == pytest.approx(4.0 * light)


class TestParamTransform:
    def test_continuous_passthrough(self):
        opt = make_optimizer([Knob("a", (0.0, 1.0)), Knob("b", (5.0, 9.0))])
        assert opt._transform_params([0.3, 7.0]) == {"a": 0.3, "b": 7.0}

    def test_discrete_rounds_and_clamps(self):
        opt = make_optimizer([Knob("m", [1.0, 2.0, 4.0, 8.0])])
        assert opt._transform_params([1.4])["m"] == 2.0
        assert opt._transform_params([99.0])["m"] == 8.0
        assert opt._transform_params([-3.0])["m"] == 1.0

    def test_bounds_for_discrete_are_index_space(self):
        opt = make_optimizer([Knob("m", [1.0, 2.0, 4.0, 8.0]), Knob("a", (2.0, 3.0))])
        assert opt._get_bounds() == [(0, 3), (2.0, 3.0)]


class TestResolveKnobs:
    def test_roles_come_from_circuit_bounds_from_config(self):
        class Toy(CircuitModel):
            KNOBS = [
                Knob("M1_GMID", role="op", sets_width_of="M1"),
                Knob("M1_L", role="geom"),
            ]

        merged = _resolve_knobs(Toy(), [Knob("M1_GMID", (5, 20)), Knob("M1_L", (1e-7, 1e-6))])
        by_name = {k.name: k for k in merged}
        assert by_name["M1_GMID"].role == "op"
        assert by_name["M1_GMID"].sets_width_of == "M1"
        assert by_name["M1_GMID"].bound == (5, 20)
        assert by_name["M1_L"].role == "geom"

    def test_missing_bound_raises(self):
        class Toy(CircuitModel):
            KNOBS = [Knob("x", role="geom")]

        with pytest.raises(ValueError, match="no bound"):
            _resolve_knobs(Toy(), [Knob("x")])


class TestWorstCase:
    def _results(self):
        return [
            CornerResult("tt", {"GBW": 12e6, "I": 8e-6, "V": 0.61, "extra": 1.0}),
            CornerResult("ss", {"GBW": 9e6, "I": 9e-6, "V": 0.55}),
            CornerResult("ff", {"GBW": 15e6, "I": 11e-6, "V": 0.64}),
        ]

    def test_max_takes_minimum(self):
        worst, binding = worst_case(self._results(), {"GBW": Spec(10e6, "max", 1.0)})
        assert worst["GBW"] == 9e6
        assert binding["GBW"] == "ss"

    def test_min_takes_maximum(self):
        worst, binding = worst_case(self._results(), {"I": Spec(10e-6, "min", 1.0)})
        assert worst["I"] == 11e-6
        assert binding["I"] == "ff"

    def test_eq_takes_furthest(self):
        worst, binding = worst_case(self._results(), {"V": Spec(0.6, "eq", 1.0)})
        assert worst["V"] == 0.55
        assert binding["V"] == "ss"

    def test_untargeted_specs_pass_through_reference(self):
        worst, _ = worst_case(self._results(), {"GBW": Spec(10e6, "max", 1.0)})
        assert worst["extra"] == 1.0

    def test_none_corners_skipped(self):
        results = [None] + self._results()
        worst, _ = worst_case(results, {"GBW": Spec(10e6, "max", 1.0)})
        assert worst["GBW"] == 9e6

    def test_empty(self):
        worst, binding = worst_case([None, None], {"GBW": Spec(1, "max", 1)})
        assert worst == {} and binding == {}


class TestComputeCostBeyondZero:
    """A 'max' spec whose actual value reaches zero or goes negative must stay graded."""

    def _opt(self, specs):
        return make_optimizer([Knob("x", (0.0, 1.0))], specs)

    def test_max_spec_graded_through_zero(self):
        opt = self._opt({"PM": Spec(60.0, "max", 1.0)})
        costs = [opt.compute_cost({"PM": a}) for a in (30.0, 10.0, 0.0, -20.0)]
        assert costs == sorted(costs)            # worse actual, higher cost, all the way down
        assert np.all(np.isfinite(costs))

    def test_max_spec_continuous_at_the_knee(self):
        opt = self._opt({"PM": Spec(60.0, "max", 1.0)})
        knee = 60.0 / np.e
        assert opt.compute_cost({"PM": knee * (1 + 1e-9)}) == pytest.approx(
            opt.compute_cost({"PM": knee * (1 - 1e-9)}), rel=1e-6)

    def test_nonpositive_target_needs_scale(self):
        with pytest.raises(ValueError, match="scale"):
            self._opt({"margin": Spec(0.0, "max", 1.0)})

    def test_scaled_spec_is_linear_in_its_unit(self):
        opt = self._opt({"margin": Spec(0.0, "max", 1.0, scale=0.05)})
        assert opt.compute_cost({"margin": 0.02}) < 1e-3
        assert opt.compute_cost({"margin": -0.10}) > opt.compute_cost({"margin": -0.05}) > 0
