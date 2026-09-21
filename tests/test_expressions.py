import numpy as np
import pytest

from mosplot.expressions import (
    Expression,
    build_expressions,
    compute_device_width,
    evaluate_expression,
)


class TestComputeDeviceWidth:
    def test_w_takes_priority(self):
        assert compute_device_width({"w": 5e-6, "weff": 1e-6}) == 5e-6

    def test_weff_scaled_by_nf(self):
        assert compute_device_width({"weff": 2e-6, "nf": 4}) == pytest.approx(8e-6)

    def test_weff_default_nf(self):
        assert compute_device_width({"weff": 2e-6}) == pytest.approx(2e-6)

    def test_missing_raises(self):
        with pytest.raises(ValueError, match="width"):
            compute_device_width({"l": 1e-6})


class TestExpression:
    def test_default_function_is_identity_for_one_variable(self):
        e = Expression(variables=["id"])
        assert e.function(3.0) == 3.0

    def test_repr_contains_variables(self):
        e = Expression(variables=["gm", "id"], label="x")
        assert "gm" in repr(e)


class TestBuildExpressions:
    @pytest.fixture(scope="class")
    def exprs(self):
        return build_expressions(width=10e-6, vdsat_var="vdsat")

    def test_gmid(self, exprs):
        e = exprs["gmid_expression"]
        assert e.variables == ["gm", "id"]
        assert e.function(20e-6, 1e-6) == pytest.approx(20.0)

    def test_current_density_captures_width(self, exprs):
        e = exprs["current_density_expression"]
        assert e.function(1e-6) == pytest.approx(0.1)

    def test_vstar(self, exprs):
        e = exprs["vstar_expression"]
        assert e.function(100e-6, 10e-6) == pytest.approx(0.2)

    def test_sign_flips(self, exprs):
        assert exprs["vsg_expression"].function(0.7) == -0.7
        assert exprs["vsb_expression"].function(-0.3) == 0.3
        assert exprs["vsd_expression"].function(0.5) == -0.5

    def test_transit_frequency(self, exprs):
        e = exprs["transit_frequency_expression"]
        assert e.function(1e-3, 1e-15) == pytest.approx(1e-3 / (2 * np.pi * 1e-15))

    def test_vdsat_var_is_respected(self):
        exprs = build_expressions(width=1e-6, vdsat_var="vdssat")
        assert exprs["vdsat_expression"].variables == ["vdssat"]


class TestEvaluateExpression:
    def test_simple_lookup(self):
        table = {"id": np.array([[1.0, 2.0], [3.0, 4.0]])}
        result, label = evaluate_expression(Expression(variables=["id"], label="L"), table)
        np.testing.assert_array_equal(result, table["id"])
        assert label == "L"

    def test_row_filter_applies_to_2d_data_only(self):
        table = {
            "gm": np.array([[1.0, 2.0], [3.0, 4.0]]),
            "id": np.array([[2.0, 2.0], [6.0, 8.0]]),
        }
        expr = Expression(variables=["gm", "id"], function=lambda g, i: g / i)
        result, _ = evaluate_expression(expr, table, filter_by_rows=np.array([1]))
        np.testing.assert_allclose(result, [[0.5, 0.5]])
