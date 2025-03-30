from typing import Any, List, Dict, Optional
from scipy.optimize import differential_evolution, OptimizeResult
from .datatypes import Spec, OptimizationParameter


class Optimizer:
    """
    Optimizes a circuit design using differential evolution.

    Attributes:
        circuit: The circuit object with an evaluate_specs method.
        parameters: A list of OptimizationParameter objects.
        target_specs: A dictionary mapping spec names to Spec objects.
        opt_params: The optimal parameters found after optimization.
        result: The optimization result from differential evolution.
    """

    def __init__(
        self,
        circuit: Any,
        parameters: List[OptimizationParameter],
        target_specs: Dict[str, Spec],
    ) -> None:
        self.circuit = circuit
        self.parameters = parameters
        self.target_specs = target_specs
        self.opt_params: Optional[Dict[str, float]] = None
        self.result: Optional[OptimizeResult] = None

    @staticmethod
    def compute_cost(specs: Dict[str, float], target_specs: Dict[str, Spec]) -> float:
        """
        Computes the cost based on actual circuit specifications and target specifications.

        Args:
            specs: A dictionary of actual circuit specifications.
            target_specs: A dictionary mapping specification names to Spec objects.

        Returns:
            The computed cost as a float.
        """
        cost = 0.0
        for key, spec_target in target_specs.items():
            target = spec_target.target
            mode = spec_target.mode
            weight = spec_target.weight
            actual = specs[key]
            error = 0.0
            if mode == "min":
                error = max(0, (actual - target) / target)
            elif mode == "max" and actual < target:
                error = max(0, (target - actual) / target)
            cost += weight * (error**2)
        return cost

    def _objective(self, x: List[float]) -> float:
        """
        Objective function for optimization.

        Args:
            x: A list of parameter values.

        Returns:
            The cost computed for the given parameter values.
        """
        params = {param.name: x[i] for i, param in enumerate(self.parameters)}
        specs = self.circuit.evaluate_specs(**params)
        return self.compute_cost(specs, self.target_specs)

    def optimize(self, maxiter: int = 10) -> OptimizeResult:
        """
        Runs the differential evolution optimization.

        Args:
            maxiter: Maximum number of iterations for the optimizer (default: 10).

        Returns:
            The optimization result.
        """
        bounds = [param.bound for param in self.parameters]
        self.result = differential_evolution(self._objective, bounds, disp=True, maxiter=maxiter, polish=True)
        self.opt_params = {
            param.name: self.result.x[i] for i, param in enumerate(self.parameters)
        }
        return self.result

    def get_opt_params(self) -> Dict[str, float]:
        """
        Returns the optimal parameters found after optimization.

        Returns:
            A dictionary of optimal parameters.
        """
        if self.opt_params is None:
            raise ValueError("Optimization has not been performed yet.")
        return self.opt_params
