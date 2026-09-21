"""FastMosfet: expression-based forward and reverse lookups."""

import numpy as np
import pytest

from mosplot.optimizer import FastMosfet

import synthetic


@pytest.fixture(scope="module")
def nmos(lookup_table):
    return FastMosfet(lookup_table, "nch", n_gmid=32, gmid_bounds=(4.0, 25.0))


class TestConstruction:
    def test_width_from_device_parameters(self, nmos):
        assert nmos._width == pytest.approx(synthetic.WIDTH)

    def test_expressions_attached(self, nmos):
        assert nmos.gmid_expression.variables == ["gm", "id"]
        assert nmos.vdsat_expression.variables == ["vdsat"]

    def test_disk_cache_roundtrip(self, lookup_table, tmp_path):
        a = FastMosfet(lookup_table, "nch", n_gmid=16, cache_dir=tmp_path)
        b = FastMosfet(lookup_table, "nch", n_gmid=16, cache_dir=tmp_path)
        args = (400e-9, 10.0, 0.6, 0.0)
        assert a.interpolate(*args, a.id_expression) == pytest.approx(
            b.interpolate(*args, b.id_expression)
        )


class TestForwardInterpolate:
    def test_single_expression(self, nmos):
        truth = synthetic.analytic_from_gmid(400e-9, 10.0, 0.6, 0.0)
        got = nmos.interpolate(400e-9, 10.0, 0.6, 0.0, nmos.id_expression)
        assert got == pytest.approx(truth["id"], rel=0.03)

    def test_expression_list(self, nmos):
        got = nmos.interpolate(
            400e-9, 10.0, 0.6, 0.0, [nmos.id_expression, nmos.gm_expression]
        )
        assert isinstance(got, list) and len(got) == 2
        assert got[1] / got[0] == pytest.approx(10.0, rel=0.02)  # gm/id

    def test_derived_expression(self, nmos):
        truth = synthetic.analytic_from_gmid(400e-9, 10.0, 0.6, 0.0)
        jd = nmos.interpolate(400e-9, 10.0, 0.6, 0.0, nmos.current_density_expression)
        assert jd == pytest.approx(truth["id"] / synthetic.WIDTH, rel=0.03)


class TestReverseLookup:
    def test_lookup_gmid_for_current_density(self, nmos):
        truth = synthetic.analytic_from_gmid(400e-9, 12.0, 0.6, 0.0)
        got = nmos.lookup_gmid_for(
            length=400e-9,
            vds=0.6,
            vbs=0.0,
            expression=nmos.current_density_expression,
            target=truth["id"] / synthetic.WIDTH,
        )
        assert got == pytest.approx(12.0, rel=0.03)

    def test_fast_method_agrees_with_scan(self, nmos):
        kwargs = dict(
            length=400e-9,
            vds=0.6,
            vbs=0.0,
            expression=nmos.current_density_expression,
            target=2.0,
        )
        scan = nmos.lookup_gmid_for(**kwargs)
        fast = nmos.lookup_gmid_for(method="fast", **kwargs)
        assert fast == pytest.approx(scan, rel=1e-6)

    def test_lookup_vds_for_id_roundtrip(self, nmos):
        # id grows (slowly) with vds through channel-length modulation, so the
        # crossing is shallow; test round-trip consistency through the same
        # table rather than against the analytic model.
        target = nmos.interpolate(400e-9, 10.0, 0.6, 0.0, nmos.id_expression)
        got = nmos.lookup_vds_for(
            length=400e-9,
            gmid=10.0,
            vbs=0.0,
            expression=nmos.id_expression,
            target=target,
        )
        assert got == pytest.approx(0.6, abs=1e-3)

    def test_axis_value_expression(self, nmos):
        # An expression of the solved axis itself: gmid such that gmid == 9.
        got = nmos.lookup_axis_for(
            "gmid",
            length=400e-9,
            vds=0.6,
            vbs=0.0,
            expression=nmos.gmid_expression,
            target=9.0,
        )
        assert got == pytest.approx(9.0, rel=0.01)

    def test_mode_all_returns_list(self, nmos):
        got = nmos.lookup_gmid_for(
            length=400e-9,
            vds=0.6,
            vbs=0.0,
            expression=nmos.current_density_expression,
            target=2.0,
            mode="all",
        )
        assert isinstance(got, list)


class TestReverseLookupValidation:
    def test_list_expression_rejected(self, nmos):
        with pytest.raises(TypeError, match="single Expression"):
            nmos.lookup_axis_for(
                "gmid",
                length=400e-9,
                vds=0.6,
                vbs=0.0,
                expression=[nmos.id_expression],
                target=1.0,
            )

    def test_unknown_axis_rejected(self, nmos):
        with pytest.raises(ValueError, match="Unknown solve axis"):
            nmos.lookup_axis_for(
                "vth",
                length=400e-9,
                vds=0.6,
                vbs=0.0,
                expression=nmos.id_expression,
                target=1.0,
            )

    def test_solved_axis_must_be_omitted(self, nmos):
        with pytest.raises(ValueError, match="must not be provided"):
            nmos.lookup_axis_for(
                "gmid",
                length=400e-9,
                vds=0.6,
                vbs=0.0,
                gmid=5.0,
                expression=nmos.id_expression,
                target=1.0,
            )

    def test_fast_mode_restrictions(self, nmos):
        common = dict(
            length=400e-9,
            vds=0.6,
            vbs=0.0,
            expression=nmos.id_expression,
            target=1.0,
            method="fast",
        )
        with pytest.raises(ValueError, match="mode='first'"):
            nmos.lookup_gmid_for(mode="all", **common)
        with pytest.raises(ValueError, match="out_of_range"):
            nmos.lookup_gmid_for(out_of_range="raise", **common)

    def test_unknown_method(self, nmos):
        with pytest.raises(ValueError, match="method"):
            nmos.lookup_gmid_for(
                length=400e-9,
                vds=0.6,
                vbs=0.0,
                expression=nmos.id_expression,
                target=1.0,
                method="banana",
            )
