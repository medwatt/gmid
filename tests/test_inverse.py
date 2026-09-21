"""Tests for the two 1-D inverse solvers.

_inverse_1d (Python) is the full-featured reference; _inverse_1d_nb (Numba)
is the optimizer hot path restricted to first-crossing + nan/clip. The last
test class checks they agree where their feature sets overlap.
"""

import numpy as np
import pytest

from mosplot.interpolation.inverse import _inverse_1d
from mosplot.interpolation.fast.interp import _inverse_1d_nb


X = np.linspace(0.0, 10.0, 101)


class TestInverse1dMonotone:
    def test_linear_crossing(self):
        y = 2.0 * X
        assert _inverse_1d(X, y, 7.0) == pytest.approx(3.5)

    def test_exact_grid_hit(self):
        y = X.copy()
        assert _inverse_1d(X, y, 4.0) == pytest.approx(4.0)

    def test_decreasing_curve(self):
        y = 10.0 - X
        assert _inverse_1d(X, y, 2.5) == pytest.approx(7.5)


class TestInverse1dOutOfRange:
    def test_nan_policy(self):
        assert np.isnan(_inverse_1d(X, X, 99.0, out_of_range="nan"))

    def test_clip_policy_returns_closest_x(self):
        assert _inverse_1d(X, X, 99.0, out_of_range="clip") == pytest.approx(10.0)

    def test_raise_policy(self):
        with pytest.raises(ValueError, match="outside"):
            _inverse_1d(X, X, 99.0, out_of_range="raise")

    def test_all_nan_curve(self):
        y = np.full_like(X, np.nan)
        assert np.isnan(_inverse_1d(X, y, 1.0))
        with pytest.raises(ValueError, match="finite"):
            _inverse_1d(X, y, 1.0, out_of_range="raise")


class TestInverse1dMultipleCrossings:
    # Parabola peaked at x=5: y = 25 - (x-5)^2 crosses 9 at x=1 and x=9.
    Y = 25.0 - (X - 5.0) ** 2

    def test_warns_and_returns_first(self):
        with pytest.warns(UserWarning, match="2 crossings"):
            assert _inverse_1d(X, self.Y, 9.0) == pytest.approx(1.0)

    def test_last(self):
        assert _inverse_1d(X, self.Y, 9.0, mode="last", warn_multiple=False) == pytest.approx(9.0)

    def test_all(self):
        sols = _inverse_1d(X, self.Y, 9.0, mode="all")
        assert sols == pytest.approx([1.0, 9.0])

    def test_nearest_picks_crossing_near_best_sample(self):
        # The sample with y closest to 24.99 is near the peak; both crossings
        # are equidistant-ish, just check a valid crossing is returned.
        sol = _inverse_1d(X, self.Y, 9.0, mode="nearest", warn_multiple=False)
        assert sol == pytest.approx(1.0) or sol == pytest.approx(9.0)

    def test_warning_suppressed(self):
        import warnings

        with warnings.catch_warnings():
            warnings.simplefilter("error")
            _inverse_1d(X, self.Y, 9.0, warn_multiple=False)


class TestInverse1dNanGaps:
    def test_crossing_not_bridged_across_gap(self):
        # y jumps across the NaN gap; without the gap there would be a
        # crossing of 5.0 between x=4 and x=6.
        y = np.where(X < 4.0, 0.0, 10.0)
        y[(X >= 4.0) & (X < 6.0)] = np.nan
        assert np.isnan(_inverse_1d(X[1:], y[1:], 5.0))

    def test_crossing_found_within_segment(self):
        y = 2.0 * X
        y[:10] = np.nan
        assert _inverse_1d(X, y, 7.0) == pytest.approx(3.5)


class TestInverse1dValidation:
    def test_bad_mode(self):
        with pytest.raises(ValueError, match="mode"):
            _inverse_1d(X, X, 1.0, mode="banana")

    def test_bad_out_of_range(self):
        with pytest.raises(ValueError, match="out_of_range"):
            _inverse_1d(X, X, 1.0, out_of_range="banana")

    def test_shape_mismatch(self):
        with pytest.raises(ValueError, match="same shape"):
            _inverse_1d(X, X[:-1], 1.0)


class TestNumbaParity:
    """_inverse_1d_nb must agree with the Python solver on its supported subset."""

    def test_linear(self):
        y = 2.0 * X
        assert _inverse_1d_nb(X, y, 7.0, 0) == pytest.approx(_inverse_1d(X, y, 7.0))

    def test_decreasing(self):
        y = 10.0 - X
        assert _inverse_1d_nb(X, y, 2.5, 0) == pytest.approx(7.5)

    def test_out_of_range_nan(self):
        assert np.isnan(_inverse_1d_nb(X, X, 99.0, 0))

    def test_out_of_range_clip(self):
        assert _inverse_1d_nb(X, X, 99.0, 1) == pytest.approx(10.0)

    def test_nan_gap_not_bridged(self):
        y = np.where(X < 4.0, 0.0, 10.0)
        y[(X >= 4.0) & (X < 6.0)] = np.nan
        assert np.isnan(_inverse_1d_nb(X[1:], y[1:], 5.0, 0))

    def test_first_crossing_of_parabola(self):
        y = 25.0 - (X - 5.0) ** 2
        assert _inverse_1d_nb(X, y, 9.0, 0) == pytest.approx(1.0)
