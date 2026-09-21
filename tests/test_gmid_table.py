"""GmIdTable: grid construction, forward/vector/curve lookups, disk cache."""

import numpy as np
import pytest

from mosplot.interpolation import GmIdTable, prebuild_fast_tables

import synthetic


GMID_BOUNDS = (4.0, 25.0)


@pytest.fixture(scope="module")
def ntable(lookup_table):
    return GmIdTable(lookup_table, "nch", n_gmid=32, gmid_bounds=GMID_BOUNDS)


class TestConstruction:
    def test_axes_are_strictly_increasing(self, ntable):
        for axis in (ntable.lengths, ntable.vbs, ntable.vds, ntable.gmid):
            assert np.all(np.diff(axis) > 0)

    def test_gmid_grid_respects_bounds(self, ntable):
        assert ntable.gmid[0] == pytest.approx(GMID_BOUNDS[0])
        assert ntable.gmid[-1] == pytest.approx(GMID_BOUNDS[1])

    def test_vgs_is_recoverable(self, ntable):
        # vgs is broadcast into the grid so a gmid query can recover it.
        assert "vgs" in ntable.params

    def test_pmos_axes_sorted_ascending(self, lookup_table):
        ptable = GmIdTable(lookup_table, "pch", n_gmid=16, gmid_bounds=GMID_BOUNDS)
        assert np.all(np.diff(ptable.vds) > 0)  # stored descending on disk
        assert ptable.vds[-1] == pytest.approx(0.0)

    def test_ambiguous_axis_lengths_raise(self, lookup_table):
        dev = dict(lookup_table["nch"])
        dev["vds"] = np.linspace(0, 1.2, len(dev["vgs"]))  # same length as vgs
        n = len(dev["vgs"])
        for p in synthetic.PARAMETER_NAMES:
            dev[p] = np.zeros((4, 3, n, n))
        with pytest.raises(ValueError, match="unidirectional"):
            GmIdTable({"bad": dev}, "bad", n_gmid=8)

    def test_invalid_gmid_bounds_raise(self, lookup_table):
        with pytest.raises(ValueError, match="bounds"):
            GmIdTable(lookup_table, "nch", n_gmid=8, gmid_bounds=(10.0, 2.0))

    def test_missing_id_or_gm_gives_clear_error(self, lookup_table):
        # A table generated without gm in parameters_to_save must fail with a
        # message naming the missing parameter, not a deep KeyError.
        dev = dict(lookup_table["nch"])
        dev.pop("gm")
        dev["parameter_names"] = [p for p in dev["parameter_names"] if p != "gm"]
        with pytest.raises(ValueError, match="parameters_to_save"):
            GmIdTable({"nch": dev}, "nch", n_gmid=8)


class TestScalarLookup:
    @pytest.mark.parametrize("gmid", [6.0, 10.0, 16.0])
    @pytest.mark.parametrize("vds", [0.3, 0.6])
    def test_id_matches_analytic(self, ntable, gmid, vds):
        truth = synthetic.analytic_from_gmid(400e-9, gmid, vds, 0.0)
        got = ntable.lookup_scalar(400e-9, gmid, vds, 0.0, params=["id"])
        assert got["id"] == pytest.approx(truth["id"], rel=0.03)

    def test_vgs_matches_analytic(self, ntable):
        truth = synthetic.analytic_from_gmid(400e-9, 10.0, 0.6, 0.0)
        got = ntable.lookup_scalar(400e-9, 10.0, 0.6, 0.0, params=["vgs"])
        assert got["vgs"] == pytest.approx(truth["vov"] + truth["vth"], rel=0.02)

    def test_body_effect(self, ntable):
        # More negative vbs raises vth, so the same gmid needs a higher vgs.
        v0 = ntable.lookup_scalar(400e-9, 10.0, 0.6, 0.0, params=["vgs"])["vgs"]
        vb = ntable.lookup_scalar(400e-9, 10.0, 0.6, -0.6, params=["vgs"])["vgs"]
        truth_shift = synthetic.vth_of_vbs(-0.6) - synthetic.vth_of_vbs(0.0)
        assert vb - v0 == pytest.approx(truth_shift, rel=0.05)

    def test_length_interpolation_between_grid_points(self, ntable):
        # 300 nm is between the 200 nm and 400 nm grid points.
        got = ntable.lookup_scalar(300e-9, 10.0, 0.6, 0.0, params=["id"])["id"]
        lo = synthetic.analytic_from_gmid(200e-9, 10.0, 0.6, 0.0)["id"]
        hi = synthetic.analytic_from_gmid(400e-9, 10.0, 0.6, 0.0)["id"]
        assert min(lo, hi) < got < max(lo, hi)


class TestVectorLookup:
    def test_matches_scalar_path(self, ntable):
        gmids = np.array([5.0, 8.0, 12.0, 20.0])
        vec = ntable.lookup(400e-9, gmids, 0.6, 0.0, params=["id", "gm"])
        for i, g in enumerate(gmids):
            scal = ntable.lookup_scalar(400e-9, g, 0.6, 0.0, params=["id", "gm"])
            assert vec["id"][i] == pytest.approx(scal["id"], rel=1e-12)
            assert vec["gm"][i] == pytest.approx(scal["gm"], rel=1e-12)

    def test_broadcasting_shape(self, ntable):
        out = ntable.lookup(400e-9, np.linspace(5, 20, 7), 0.6, 0.0, params=["id"])
        assert out["id"].shape == (7,)


class TestCurveLookup:
    def test_curve_matches_scalar_lookups(self, ntable):
        axis, curves = ntable.lookup_curve(
            solve_for="gmid", params=["id"], length=400e-9, vds=0.6, vbs=0.0
        )
        np.testing.assert_array_equal(axis, ntable.gmid)
        for i in [0, 10, 20, 31]:
            scal = ntable.lookup_scalar(400e-9, axis[i], 0.6, 0.0, params=["id"])
            assert curves["id"][i] == pytest.approx(scal["id"], rel=1e-9)

    def test_solved_axis_must_not_be_given(self, ntable):
        with pytest.raises(ValueError, match="must not be provided"):
            ntable.lookup_curve(
                solve_for="gmid", params=["id"], gmid=10.0, length=400e-9, vds=0.6, vbs=0.0
            )

    def test_missing_fixed_axis_raises(self, ntable):
        with pytest.raises(ValueError, match="Missing fixed axes"):
            ntable.lookup_curve(solve_for="gmid", params=["id"], length=400e-9, vds=0.6)

    def test_unknown_param_raises(self, ntable):
        with pytest.raises(KeyError, match="banana"):
            ntable.lookup_curve(
                solve_for="gmid", params=["banana"], length=400e-9, vds=0.6, vbs=0.0
            )


class TestDiskCache:
    def test_save_load_roundtrip(self, ntable, tmp_path):
        path = tmp_path / "cache.npz"
        ntable.save(path)
        loaded = GmIdTable.load(path)
        assert all(np.shares_memory(a, loaded._data_stack) for a in loaded._data.values())
        a = ntable.lookup_scalar(400e-9, 10.0, 0.6, 0.0)
        b = loaded.lookup_scalar(400e-9, 10.0, 0.6, 0.0)
        for k in a:
            assert b[k] == pytest.approx(a[k], rel=1e-12, abs=1e-300, nan_ok=True)

    def test_build_or_load_creates_then_reuses(self, lookup_table, tmp_path, capsys):
        kwargs = dict(cache_dir=tmp_path, n_gmid=16, gmid_bounds=GMID_BOUNDS, verbose=True)
        GmIdTable.build_or_load(lookup_table, "nch", **kwargs)
        assert "Building" in capsys.readouterr().out
        GmIdTable.build_or_load(lookup_table, "nch", **kwargs)
        assert "Loading cached" in capsys.readouterr().out

    def test_cache_invalidated_when_data_changes(self, lookup_table, tmp_path):
        kwargs = dict(cache_dir=tmp_path, n_gmid=16, gmid_bounds=GMID_BOUNDS)
        GmIdTable.build_or_load(lookup_table, "nch", **kwargs)
        n_before = len(list(tmp_path.glob("*.npz")))

        modified = {"nch": dict(lookup_table["nch"])}
        modified["nch"]["id"] = lookup_table["nch"]["id"] * 1.1
        GmIdTable.build_or_load(modified, "nch", **kwargs)
        # The data fingerprint changed, so a second cache file must appear.
        assert len(list(tmp_path.glob("*.npz"))) == n_before + 1

    def test_prebuild_returns_existing_paths(self, lookup_table, tmp_path):
        paths = prebuild_fast_tables(
            lookup_table,
            cache_dir=tmp_path,
            device_keys=["nch"],
            n_gmid=16,
            gmid_bounds=GMID_BOUNDS,
            verbose=False,
        )
        assert paths["nch"].exists()

    def test_prebuild_accepts_npz_path(self, lut_npz_path, tmp_path):
        paths = prebuild_fast_tables(
            str(lut_npz_path),
            cache_dir=tmp_path,
            device_keys=["pch"],
            n_gmid=16,
            gmid_bounds=GMID_BOUNDS,
            verbose=False,
        )
        assert paths["pch"].exists()
