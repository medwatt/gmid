"""DC solver + sizing/recornering flow on an analytic toy model.

The toy "circuit" has a closed-form solution, so the test can check that:
  * solve_point drives the residuals to zero at the known solution,
  * size_reference freezes the correct geometry,
  * recorner re-solves the op knob so the frozen width is reproduced on a
    corner with shifted "process" (the cond["k"] factor),
  * RECORNER_RESOLVE knobs are re-solved to conserve a derived quantity.

Device model: jd = k * G**2 (k plays the role of process), W = IB / jd,
and one unknown VX pinned by the residual VX == 2 / G.
"""

from types import SimpleNamespace

import pytest

from mosplot.optimizer import CircuitModel, Knob, State, Unknown, vres, rres
from mosplot.optimizer.evaluate import check_square, evaluate_corners, recorner, size_reference
from mosplot.optimizer.solve import solve_point


def corner(name, k, iref=10e-6):
    return SimpleNamespace(name=name, conditions={"k": k, "iref": iref})


class Toy(CircuitModel):
    KNOBS = [
        Knob("G", (2.0, 30.0), role="op", sets_width_of="M1"),
        Knob("L1", (1e-7, 1e-6), role="geom"),
        Knob("IB", (1e-6, 1e-4), role="external"),
    ]
    UNKNOWNS = [Unknown("VX", seed=lambda c: 1.0, bound=(0.01, 5.0))]

    def solve_point(self, v, dev, cond):
        jd = cond["k"] * v.G**2
        W = {"M1": v.IB / jd}
        L = {"M1": v.L1}
        return State(W=W, L=L, G=v.G, VX=v.VX, IB=v.IB, k=cond["k"])

    def residuals(self, b):
        return [vres(b.VX, 2.0 / b.G)]

    def specs(self, b, cond):
        return {"W": b.W["M1"], "VX": b.VX, "G": b.G}


class ToyWithConservation(Toy):
    """Adds a derived conserved quantity: IREF = IB * k, re-solved via IB."""

    RECORNER_RESOLVE = ["IB"]

    def solve_point(self, v, dev, cond):
        b = super().solve_point(v, dev, cond)
        b.IREF = v.IB * cond["k"]
        return b

    def freeze_extra(self, b):
        return {"IREF": b.IREF}

    def recorner_residuals(self, b, frozen):
        e = frozen["extra"]
        return [rres(b.IREF, e["IREF"], e["IREF"])]


KNOBS = {"G": 10.0, "L1": 4e-7, "IB": 10e-6}


class TestSolvePoint:
    def test_unknown_converges_to_closed_form(self):
        b, v, r = solve_point(Toy(), KNOBS, corner("tt", k=1.0))
        assert b.VX == pytest.approx(0.2, abs=1e-9)
        assert r < 1e-9

    def test_no_free_names_skips_solver(self):
        b, v, r = solve_point(Toy(), dict(KNOBS, VX=0.5), corner("tt", 1.0), free_names=[])
        assert b.VX == 0.5 and r == 0.0

    def test_pinned_overrides(self):
        b, _, _ = solve_point(
            Toy(), KNOBS, corner("tt", 1.0), free_names=[], pinned={"VX": 0.33}
        )
        assert b.VX == 0.33

    def test_seed_is_clipped_into_bounds(self):
        b, _, r = solve_point(Toy(), KNOBS, corner("tt", 1.0), seed={"VX": 100.0})
        assert b.VX == pytest.approx(0.2, abs=1e-6)


class TestSizeReference:
    def test_freezes_geometry_and_knobs(self):
        res, frozen = size_reference(Toy(), KNOBS, corner("tt", k=1.0))
        assert res.specs["W"] == pytest.approx(1e-7)  # IB / (k G^2) = 1e-5/100
        assert frozen["W"]["M1"] == pytest.approx(1e-7)
        assert frozen["L"]["M1"] == 4e-7
        assert frozen["knobs"] == KNOBS
        assert frozen["unknown_seed"]["VX"] == pytest.approx(0.2, abs=1e-6)
        assert frozen["dims"]["M1"]["Width"] == pytest.approx(1e-7)


class TestRecorner:
    def test_op_knob_resolved_to_match_width(self):
        _, frozen = size_reference(Toy(), KNOBS, corner("tt", k=1.0))
        res = recorner(Toy(), frozen, corner("ss", k=1.21))
        # W = IB/(k G^2) conserved with IB pinned -> G scales by 1/sqrt(k).
        assert res.specs["G"] == pytest.approx(10.0 / 1.1, rel=1e-6)
        assert res.specs["W"] == pytest.approx(frozen["W"]["M1"], rel=1e-9)
        # The unknown tracks the re-solved knob.
        assert res.specs["VX"] == pytest.approx(2.0 / (10.0 / 1.1), rel=1e-6)

    def test_recorner_resolve_conserves_derived_quantity(self):
        model = ToyWithConservation()
        _, frozen = size_reference(model, KNOBS, corner("tt", k=1.0))
        res = recorner(model, frozen, corner("ss", k=1.21))
        # IREF = IB*k conserved -> IB drops by 1/1.21; W = IB/(k G^2) frozen
        # -> k*G^2 must drop by 1.21 -> G = 10/1.21.
        assert res.op.IB == pytest.approx(10e-6 / 1.21, rel=1e-6)
        assert res.op.IREF == pytest.approx(10e-6, rel=1e-9)
        assert res.specs["G"] == pytest.approx(10.0 / 1.21, rel=1e-6)
        assert res.specs["W"] == pytest.approx(frozen["W"]["M1"], rel=1e-9)


class TestEvaluateCorners:
    def test_reference_plus_recorners(self):
        corners = [corner("tt", 1.0), corner("ss", 1.21), corner("ff", 0.81)]
        results = evaluate_corners(Toy(), KNOBS, corners)
        assert [r.name for r in results] == ["tt", "ss", "ff"]
        for r in results:
            assert r.specs["W"] == pytest.approx(1e-7, rel=1e-8)
        assert results[1].specs["G"] == pytest.approx(10.0 / 1.1, rel=1e-6)
        assert results[2].specs["G"] == pytest.approx(10.0 / 0.9, rel=1e-6)

    def test_nonzero_ref_index(self):
        corners = [corner("ss", 1.21), corner("tt", 1.0)]
        results = evaluate_corners(Toy(), KNOBS, corners, ref_index=1)
        assert results[1].specs["G"] == 10.0  # reference keeps the knob value

    def test_failed_reference_returns_all_none(self):
        class Broken(Toy):
            def solve_point(self, v, dev, cond):
                raise RuntimeError("boom")

        results = evaluate_corners(Broken(), KNOBS, [corner("tt", 1.0), corner("ss", 1.2)])
        assert results == [None, None]

    def test_failed_single_corner_is_none(self):
        class BreaksOffReference(Toy):
            def solve_point(self, v, dev, cond):
                if cond["k"] != 1.0:
                    raise RuntimeError("boom")
                return super().solve_point(v, dev, cond)

        results = evaluate_corners(
            BreaksOffReference(), KNOBS, [corner("tt", 1.0), corner("ss", 1.2)]
        )
        assert results[0] is not None
        assert results[1] is None


class TestRecornerBounds:
    def test_op_knob_may_leave_its_search_box_on_another_corner(self):
        # G searched in (2, 30); sized at 28 in tt, the frozen width needs G = 28/0.9 in ff.
        knobs = dict(KNOBS, G=28.0)
        _, frozen = size_reference(Toy(), knobs, corner("tt", k=1.0))
        res = recorner(Toy(), frozen, corner("ff", k=0.81))
        assert res.max_residual < 1e-9
        assert res.specs["G"] == pytest.approx(28.0 / 0.9, rel=1e-6)


class TestSquareCheck:
    def test_conserved_quantity_without_resolve_knob_is_rejected(self):
        class Unbalanced(ToyWithConservation):
            RECORNER_RESOLVE = []

        corners = [corner("tt", 1.0), corner("ss", 1.21)]
        with pytest.raises(ValueError, match="RECORNER_RESOLVE"):
            check_square(Unbalanced(), KNOBS, corners)
        assert check_square(ToyWithConservation(), KNOBS, corners)["recorner"] == (3, 3)
