from scipy.optimize import differential_evolution

class Optimizer:
    def __init__(self, circuit, var_names, bounds, target_specs):
        self.circuit = circuit
        self.var_names = var_names
        self.bounds = bounds
        self.target_specs = target_specs
        self.opt_params = None
        self.result = None

    @staticmethod
    def compute_cost(specs, target_specs):
        cost = 0.0
        for key, t in target_specs.items():
            target = t["target"]
            mode = t["mode"]
            weight = t["weight"]
            actual = specs[key]
            error = 0
            if mode == "min":
                error = max(0, (actual - target) / target)
            elif mode == "max" and actual < target:
                error = max(0, (target - actual) / target)
            cost += weight * (error ** 2)
        return cost

    def _objective(self, x):
        params = {name: x[i] for i, name in enumerate(self.var_names)}
        specs = self.circuit.evaluate_specs(**params)
        return self.compute_cost(specs, self.target_specs)

    def optimize(self, maxiter=10):
        self.result = differential_evolution(self._objective, self.bounds, disp=True, maxiter=maxiter, polish=True)
        self.opt_params = {name: self.result.x[i] for i, name in enumerate(self.var_names)}
        return self.result

    def get_opt_params(self):
        return self.opt_params

    def get_opt_specs(self):
        return self.circuit.evaluate_specs(**self.opt_params)
