from __future__ import annotations

from typing import Any
from operator import index
import numpy as np
import numpy.typing as npt
import cma
from scipy.optimize import OptimizeResult, minimize
from .evaluate import (
    CornerResult,
    check_square,
    evaluate_corners,
    recorner,
    size_reference,
    worst_case,
)
from .datatypes import Knob, Spec


_SOFTPLUS_LIN_THRESHOLD = 20.0


def _softplus(x: float, beta: float = 20.0) -> float:
    """Smooth, differentiable approximation of max(0, x)."""
    if x > _SOFTPLUS_LIN_THRESHOLD / beta:
        return x
    return float(np.log1p(np.exp(beta * x)) / beta)


def _resolve_knobs(model, config_knobs: list[Knob]) -> list[Knob]:
    """Merge the circuit's structural knobs (roles, no bounds) with the config's bounds.

    Roles and sets_width_of come from the circuit (physics); bounds come from the config
    (PDK/spec). This is what lets a circuit file carry zero numeric constants.
    """
    roles = {k.name: (k.role, k.sets_width_of, k.recorner_bound) for k in model.KNOBS}
    merged = []
    for ck in config_knobs:
        if ck.bound is None:
            raise ValueError(f"knob {ck.name!r} has no bound; set it in the config PARAMETERS.")
        role, sets_width_of, rb = roles.get(ck.name, (ck.role, ck.sets_width_of, None))
        merged.append(Knob(ck.name, ck.bound, role, sets_width_of, ck.recorner_bound or rb))
    return merged


class Optimizer:
    """Optimizes a circuit design using CMA-ES followed by SLSQP polishing."""

    def __init__(
        self,
        model: Any,
        parameters: list[Knob] | None,
        target_specs: dict[str, Spec],
        corners: list[Any],
        ref_index: int = 0,
        executor: Any = None,
    ) -> None:
        self.model = model
        # Merge config bounds onto the circuit's structural roles, then make the resolved
        # knobs the single source of truth for both the optimizer and the evaluator.
        self.parameters = _resolve_knobs(
            model, parameters if parameters is not None else model.KNOBS
        )
        model.KNOBS = self.parameters
        for key, spec in target_specs.items():
            if spec.scale is None and (spec.target <= 0 if spec.mode != "eq" else spec.target == 0):
                raise ValueError(f"spec {key!r}: a target of {spec.target} has no relative scale; "
                                 f"give the Spec a `scale` (in the spec's unit).")
        self.target_specs = target_specs
        self.corners = corners
        self.ref_index = ref_index
        # A non-square multicorner re-solve fails silently (every corner looks infeasible):
        # check the equation count once, at the middle of the search box (of a discrete set).
        mid = {k.name: (0.5 * (k.bound[0] + k.bound[1]) if isinstance(k.bound, tuple)
                        else float(np.asarray(k.bound)[len(k.bound) // 2])) for k in self.parameters}
        check_square(model, mid, corners, ref_index)
        self.executor = executor
        self.opt_params: dict[str, float] | None = None
        self.result: OptimizeResult | None = None
        self.corner_results: list[CornerResult] | None = None
        self.binding: dict[str, str] | None = None
        # Populated by _run_corner_analysis for the netlist writer / report.
        self.frozen: dict | None = None
        self.reference_op: Any = None

    def compute_cost(self, specs: dict[str, float]) -> float:
        """
        Scalar cost for a given set of circuit specs against the targets.

        "max" specs use log-scale normalisation so that large violations
        (e.g. CMRR 10x below target) produce proportionally larger gradients
        than a linearisation would. "min" specs use linear normalisation,
        which is appropriate for quantities that vary by small factors, and
        include a secondary term that keeps cost non-zero when satisfied so
        the optimizer continues driving them down. "eq" specs use squared
        relative error, penalising deviation from the target in both directions.
        """
        cost = 0.0
        for key, spec in self.target_specs.items():
            actual = specs.get(key)
            if actual is None:
                cost += 1e6
                continue
            target = spec.target

            relative = spec.scale is None      # validated in __init__: then target > 0 or eq != 0
            scale = abs(target) if relative else spec.scale

            if spec.mode == "min":
                raw = (actual - target) / scale
                violation = _softplus(raw)
                secondary = max(0.0, 1.0 + raw)      # = actual/target for a relative spec
                cost += spec.weight * (violation**2 + 0.05 * secondary)

            elif spec.mode == "max":
                if relative:
                    # log(target/actual), continued linearly (value and slope matched) below
                    # target/e, so a zero or negative actual is graded instead of a cliff
                    knee = target / np.e
                    raw = (float(np.log(target / actual)) if actual >= knee
                           else 2.0 - np.e * actual / target)
                else:
                    raw = (target - actual) / scale
                violation = _softplus(raw)
                cost += spec.weight * (violation**2)

            elif spec.mode == "eq":
                raw = (actual - target) / scale
                cost += spec.weight * (raw**2)

        return cost

    def _transform_params(self, x: list[float] | npt.NDArray) -> dict[str, float]:
        params: dict[str, float] = {}
        for i, param in enumerate(self.parameters):
            bound = param.bound
            if isinstance(bound, (list, np.ndarray)):
                idx = max(0, min(int(round(x[i])), len(bound) - 1))
                params[param.name] = float(np.array(bound)[idx])
            else:
                params[param.name] = float(x[i])
        return params

    # Infeasible designs return a penalty ABOVE this floor, GRADED by how badly the DC solve
    # missed convergence. A flat penalty would make every failed candidate identical, leaving
    # CMA-ES unable to rank a failed population and stalling it on a plateau; grading by the
    # residual gives a gradient back toward the region where the operating point solves.
    _INFEASIBLE_FLOOR = 1e6

    def _objective(self, x: list[float] | npt.NDArray) -> float:
        knobs = self._transform_params(x)
        results = evaluate_corners(self.model, knobs, self.corners, self.ref_index)

        infeasible = 0.0
        for r in results:
            if r is None:
                infeasible += 10.0  # solve raised / no operating point
            elif r.max_residual > 1e-3:
                infeasible += min(r.max_residual, 10.0)  # solver did not converge here
        if infeasible > 0.0:
            return self._INFEASIBLE_FLOOR * (1.0 + infeasible)

        worst, _ = worst_case(results, self.target_specs)
        if not worst:
            return self._INFEASIBLE_FLOOR
        return self.compute_cost(worst)

    def _get_bounds(self) -> list[tuple[float, float]]:
        bounds = []
        for param in self.parameters:
            bound = param.bound
            if isinstance(bound, (list, np.ndarray)):
                bounds.append((0, len(bound) - 1))
            else:
                bounds.append(bound)
        return bounds

    def _run_cma(
        self,
        x0_norm: npt.NDArray,
        lbs: npt.NDArray,
        ubs: npt.NDArray,
        maxiter: int,
        sigma0: float,
        seed: int,
        pool=None,
    ) -> tuple[npt.NDArray, float, int]:
        n = len(lbs)

        def objective_normalised(x_norm: npt.NDArray) -> float:
            return self._objective(lbs + x_norm * (ubs - lbs))

        opts = cma.CMAOptions()
        opts["bounds"] = [[0.0] * n, [1.0] * n]
        opts["maxiter"] = maxiter
        opts["seed"] = seed
        opts["verb_log"] = 0
        opts["verbose"] = 1
        opts["tolx"] = 1e-5
        opts["tolfun"] = 1e-5

        es = cma.CMAEvolutionStrategy(x0_norm, sigma0, opts)
        if pool is None:
            es.optimize(objective_normalised)
        else:
            from .parallel import objective

            while not es.stop() or es.countiter == 0:
                candidates = es.ask()
                costs = list(pool.map(objective, [lbs + x * (ubs - lbs) for x in candidates]))
                es.tell(candidates, costs)
                es.disp()

        best_x = lbs + es.result.xbest * (ubs - lbs)
        return best_x, float(es.result.fbest), int(es.result.iterations)

    def _run_corner_analysis(self) -> None:
        """Re-evaluate every corner at the optimum, freezing the reference geometry, and
        store the per-corner results plus the frozen design for the netlist writer."""
        opt_params = self.get_opt_params()
        ref = self.corners[self.ref_index]
        ref_res, frozen = size_reference(self.model, opt_params, ref)

        results: list[CornerResult | None] = [None] * len(self.corners)
        results[self.ref_index] = ref_res
        for i, corner in enumerate(self.corners):
            if i == self.ref_index:
                continue
            results[i] = recorner(self.model, frozen, corner)

        for r in results:
            if r is not None:
                r.cost = self.compute_cost(r.specs) if r.specs else float("inf")

        self.corner_results = results
        self.frozen = frozen
        self.reference_op = ref_res.op
        _, self.binding = worst_case(results, self.target_specs)

    def optimize(
        self,
        maxiter: int = 500,
        sigma0: float = 0.3,
        n_restarts: int = 1,
        seed: int = 42,
        workers: int = 1,
    ) -> OptimizeResult:
        """
        Run CMA-ES (with optional random restarts) then polish with SLSQP.

        Args:
            maxiter: Maximum CMA-ES iterations per restart.
            sigma0: Initial step size as a fraction of the normalised [0, 1] search space.
            n_restarts: Number of independent CMA-ES runs. The first starts from the
                centre of the search space; subsequent ones from random points. Increase
                beyond 1 only if results are inconsistent across runs (multimodal landscape).
            seed: Base random seed; restart i uses seed + i for reproducibility.
            workers: CMA-ES worker processes. 1 runs serially; larger values share lookup
                arrays using spawn. Call under an ``if __name__ == "__main__":`` guard.
        """
        from .parallel import objective_pool

        if isinstance(workers, bool):
            raise ValueError("workers must be a positive integer")
        try:
            workers = index(workers)
        except TypeError as exc:
            raise ValueError("workers must be a positive integer") from exc
        if workers < 1:
            raise ValueError("workers must be a positive integer")
        bounds = self._get_bounds()
        lbs = np.array([b[0] for b in bounds], dtype=float)
        ubs = np.array([b[1] for b in bounds], dtype=float)
        rng = np.random.default_rng(seed)

        best_x: npt.NDArray = np.empty(len(bounds))
        best_cost = np.inf
        total_iters = 0

        with objective_pool(self, workers) as pool:
            for i in range(n_restarts):
                x0_norm = np.full(len(bounds), 0.5) if i == 0 else rng.uniform(0.1, 0.9, len(bounds))
                print(f"\n--- CMA-ES restart {i + 1}/{n_restarts} ---")
                x, cost, iters = self._run_cma(x0_norm, lbs, ubs, maxiter, sigma0, seed + i, pool)
                total_iters += iters
                if cost < best_cost:
                    best_cost = cost
                    best_x = x

        print(f"\nBest CMA-ES cost across {n_restarts} restart(s): {best_cost:.6g}")
        print("Polishing with SLSQP...")

        polish = minimize(
            self._objective,
            best_x,
            method="SLSQP",
            bounds=bounds,
            options={"maxiter": 500, "ftol": 1e-9},
        )

        if polish.fun < best_cost:
            best_x, best_cost = polish.x, polish.fun
            print(f"SLSQP improved cost to {best_cost:.6g}")
        else:
            print("SLSQP did not improve on CMA-ES result.")

        self.opt_params = self._transform_params(best_x)
        self.result = OptimizeResult(
            x=best_x,
            fun=best_cost,
            success=True,
            message="CMA-ES + SLSQP optimization complete.",
            nit=total_iters,
        )

        self._run_corner_analysis()

        return self.result

    def get_opt_params(self) -> dict[str, float]:
        if self.opt_params is None:
            raise ValueError("Optimization has not been performed yet.")
        return self.opt_params
