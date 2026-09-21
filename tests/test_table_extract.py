"""extract_2d_table / tile_arrays: slicing the raw 4-D table for plotting."""

import numpy as np
import pytest

from mosplot.table import extract_2d_table
from mosplot.table.extract import tile_arrays

import synthetic


@pytest.fixture(scope="module")
def dev(lookup_table):
    return lookup_table["nch"]


N_L = len(synthetic.LENGTHS)
N_VBS = len(synthetic.VBS_AXIS)
N_VGS = len(synthetic.VGS_AXIS)
N_VDS = len(synthetic.VDS_AXIS)


class TestTileArrays:
    def test_axis_matches_rows(self):
        a = np.arange(3.0)
        b = np.zeros((3, 5))
        ta, tb = tile_arrays(a, b)
        assert ta.shape == (3, 5)
        np.testing.assert_array_equal(ta[:, 0], a)

    def test_axis_matches_cols(self):
        a = np.arange(5.0)
        b = np.zeros((3, 5))
        ta, _ = tile_arrays(a, b)
        assert ta.shape == (3, 5)
        np.testing.assert_array_equal(ta[0], a)

    def test_2d_2d_warns(self):
        with pytest.warns(UserWarning, match="unexpected shapes"):
            tile_arrays(np.zeros((2, 2)), np.zeros((3, 3)))


class TestExtract2dTable:
    def test_requires_two_fixed_axes(self, dev):
        with pytest.raises(ValueError, match="at least two"):
            extract_2d_table(lookup_table=dev, vds=0.6)

    def test_fixed_scalars_snap_to_grid(self, dev):
        _, filt, table = extract_2d_table(
            lookup_table=dev, vbs=0.0, vds=0.61, vgs=(0.0, 1.2)
        )
        assert filt["vds"] == pytest.approx(0.6)  # nearest grid point
        assert filt["vbs"] == pytest.approx(0.0)
        # Remaining free axes: length (4) x vgs (61)
        assert table["id"].shape == (N_L, N_VGS)

    def test_secondary_variable_detection(self, dev):
        secondary, _, _ = extract_2d_table(
            lookup_table=dev, vbs=0.0, vds=0.6, vgs=(0.0, 1.2), primary="vgs"
        )
        assert secondary == "length"

    def test_range_with_step(self, dev):
        _, filt, _ = extract_2d_table(
            lookup_table=dev, vbs=0.0, vds=0.6, vgs=(0.0, 1.2, 0.1), primary="vgs"
        )
        # step 0.1 on a 0.02 grid -> every 5th point, endpoint kept
        assert filt["vgs"][0] == pytest.approx(0.0)
        assert filt["vgs"][-1] == pytest.approx(1.2)
        assert filt["vgs"][1] - filt["vgs"][0] == pytest.approx(0.1)

    def test_length_list_selection(self, dev):
        _, filt, table = extract_2d_table(
            lookup_table=dev, length=[100e-9, 400e-9], vbs=0.0, vds=0.6
        )
        np.testing.assert_allclose(filt["length"], [100e-9, 400e-9])
        assert table["id"].shape[0] in (2, N_VGS)  # 2 rows after orientation

    def test_extracted_values_match_source(self, dev):
        _, filt, table = extract_2d_table(lookup_table=dev, vbs=0.0, vds=0.6)
        # id at (L=400n, vbs=0, vgs=0.8, vds=0.6) must survive the slicing.
        li = 2  # 400 nm
        gi = int(np.argmin(np.abs(synthetic.VGS_AXIS - 0.8)))
        di = int(np.argmin(np.abs(synthetic.VDS_AXIS - 0.6)))
        expected = dev["id"][li, 2, gi, di]
        assert table["id"][li, gi] == pytest.approx(expected)

    def test_axes_are_tiled_to_data_shape(self, dev):
        _, _, table = extract_2d_table(lookup_table=dev, vbs=0.0, vds=0.6)
        assert table["vgs"].shape == table["id"].shape
        assert table["length"].shape == table["id"].shape
