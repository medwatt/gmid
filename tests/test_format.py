import pytest

from mosplot.util.format import si, si_area


class TestSi:
    def test_none(self):
        assert si(None) == "N/A"

    @pytest.mark.parametrize(
        "value, expected",
        [
            (0.0, "0"),
            (1.0, "1"),
            (1234.0, "1.234k"),
            (10e6, "10M"),
            (2.5e9, "2.5G"),
            (3e12, "3T"),
            (1e-3, "1m"),
            (47e-6, "47u"),
            (1e-9, "1n"),
            (2.2e-12, "2.2p"),
            (5e-15, "5f"),
            (-47e-6, "-47u"),
        ],
    )
    def test_ascii_prefixes(self, value, expected):
        assert si(value) == expected

    def test_unicode_micro(self):
        assert si(47e-6, unicode=True) == "47µ"

    def test_rounding_steps_up_a_tier(self):
        # 999.96 rounds to 1000 at 4 sig figs -> must promote to the next prefix.
        assert si(999.96e-6) == "1m"

    def test_tiny_value_falls_through(self):
        out = si(1e-20)
        assert "e" in out  # below femto: plain scientific notation


class TestSiArea:
    def test_none(self):
        assert si_area(None) == "N/A"

    @pytest.mark.parametrize(
        "value, expected",
        [
            (2.0, "2 m²"),
            (5e-12, "5 µm²"),
            (50e-12, "50 µm²"),
            (3e-18, "3 nm²"),
        ],
    )
    def test_units(self, value, expected):
        assert si_area(value) == expected
