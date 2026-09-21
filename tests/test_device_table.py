"""DeviceTable / DevicePoint: the polarity-normalising lookup API circuits use."""

import pytest

from mosplot.optimizer import DeviceTable, FastMosfet

import synthetic


GMID_BOUNDS = (4.0, 25.0)


@pytest.fixture(scope="module")
def nmos(lookup_table):
    return DeviceTable(FastMosfet(lookup_table, "nch", gmid_bounds=GMID_BOUNDS), "n")


@pytest.fixture(scope="module")
def pmos(lookup_table):
    return DeviceTable(FastMosfet(lookup_table, "pch", gmid_bounds=GMID_BOUNDS), "p")


class TestPolarity:
    def test_invalid_polarity(self, lookup_table):
        with pytest.raises(AssertionError):
            DeviceTable(FastMosfet(lookup_table, "nch", gmid_bounds=GMID_BOUNDS), "x")

    def test_nmos_point_matches_analytic(self, nmos):
        truth = synthetic.analytic_from_gmid(400e-9, 10.0, 0.6, 0.0)
        pt = nmos(gmid=10.0, L=400e-9, vds=0.6, vsb=0.0)
        assert pt.vgs == pytest.approx(truth["vov"] + truth["vth"], rel=0.02)
        assert pt.vdsat == pytest.approx(truth["vov"], rel=0.02)
        assert pt.jd == pytest.approx(truth["id"] / synthetic.WIDTH, rel=0.03)
        assert pt.gds_id == pytest.approx(truth["gds"] / truth["id"], rel=0.03)
        assert pt.vds_used == 0.6

    def test_pmos_point_is_positive_magnitudes(self, pmos):
        # The PMOS table stores vgs/vdsat negative; the DevicePoint API
        # exposes positive magnitudes for both polarities.
        truth = synthetic.analytic_from_gmid(400e-9, 10.0, 0.6, 0.0)
        pt = pmos(gmid=10.0, L=400e-9, vds=0.6, vsb=0.0)
        assert pt.vgs == pytest.approx(truth["vov"] + truth["vth"], rel=0.02)
        assert pt.vdsat == pytest.approx(truth["vov"], rel=0.02)
        assert pt.jd > 0.0

    def test_nmos_pmos_symmetric_table_agree(self, nmos, pmos):
        n = nmos(gmid=8.0, L=200e-9, vds=0.5, vsb=0.0)
        p = pmos(gmid=8.0, L=200e-9, vds=0.5, vsb=0.0)
        assert p.vgs == pytest.approx(n.vgs, rel=1e-6)
        assert p.jd == pytest.approx(n.jd, rel=1e-6)

    def test_body_effect_through_vsb(self, nmos):
        pt0 = nmos(gmid=10.0, L=400e-9, vds=0.6, vsb=0.0)
        ptb = nmos(gmid=10.0, L=400e-9, vds=0.6, vsb=0.6)
        shift = synthetic.vth_of_vbs(-0.6) - synthetic.vth_of_vbs(0.0)
        assert ptb.vgs - pt0.vgs == pytest.approx(shift, rel=0.05)

    def test_gmb_present_when_table_has_gmbs(self, nmos):
        pt = nmos(gmid=10.0, L=400e-9, vds=0.6)
        assert pt.gmb_id == pytest.approx(0.2 * 10.0, rel=0.02)  # gmbs = 0.2*gm


class TestSmallSignal:
    def test_scaling_to_real_width(self, nmos):
        pt = nmos(gmid=10.0, L=400e-9, vds=0.6)
        idd = 20e-6
        width = 2.0 * nmos.table_width
        ss = pt.small_signal(10.0, idd, width, nmos.table_width)
        assert ss["gm"] == pytest.approx(10.0 * idd)
        assert ss["gds"] == pytest.approx(pt.gds_id * idd)
        assert ss["cgs"] == pytest.approx(pt.cgs * 2.0)
        assert ss["cgd"] == pytest.approx(pt.cgd * 2.0)
        assert ss["cdd"] == pytest.approx(pt.cdd * 2.0)

    def test_gmb_toggle(self, nmos):
        pt = nmos(gmid=10.0, L=400e-9, vds=0.6)
        on = pt.small_signal(10.0, 1e-6, 1e-6, nmos.table_width, use_gmb=True)
        off = pt.small_signal(10.0, 1e-6, 1e-6, nmos.table_width, use_gmb=False)
        assert on["gmb"] > 0.0
        assert off["gmb"] == 0.0


class TestGmbAbsentTable:
    def test_gmb_zero_when_lut_lacks_gmbs(self, lookup_table):
        dev = dict(lookup_table["nch"])
        dev.pop("gmbs")
        dev["parameter_names"] = [p for p in dev["parameter_names"] if p != "gmbs"]
        with pytest.warns(UserWarning, match="'nch' has no gmbs"):
            table = DeviceTable(
                FastMosfet({"nch": dev}, "nch", gmid_bounds=GMID_BOUNDS), "n"
            )
        pt = table(gmid=10.0, L=400e-9, vds=0.6)
        assert pt.gmb_id == 0.0


class TestPmosBodyEffectAxisConvention:
    """The LUT generator stores the PMOS vbs axis >= 0 (bulk above source); reverse body bias
    must be looked up on that side of the axis."""

    @pytest.fixture(scope="class")
    def pmos_up(self):
        table = {"pch": synthetic.build_device_vbs_up("p"), "description": "", "simulator": "synthetic",
                 "parameter_names": list(synthetic.PARAMETER_NAMES), "device_parameters": {"w": synthetic.WIDTH}}
        return DeviceTable(FastMosfet(table, "pch", gmid_bounds=GMID_BOUNDS), "p")

    def test_body_effect_with_positive_vbs_axis(self, pmos_up):
        pt0 = pmos_up(gmid=10.0, L=400e-9, vds=0.6, vsb=0.0)
        ptb = pmos_up(gmid=10.0, L=400e-9, vds=0.6, vsb=0.6)
        shift = synthetic.vth_of_vbs(-0.6) - synthetic.vth_of_vbs(0.0)
        assert ptb.vgs - pt0.vgs == pytest.approx(shift, rel=0.05)

    def test_same_answer_as_negative_axis_table(self, pmos, pmos_up):
        a = pmos(gmid=8.0, L=200e-9, vds=0.5, vsb=0.3)
        b = pmos_up(gmid=8.0, L=200e-9, vds=0.5, vsb=0.3)
        assert b.vgs == pytest.approx(a.vgs, rel=1e-6)
