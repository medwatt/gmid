# imports <<<
import os

import numpy as np

from .simulators.spice_simulators.mosfet_simulation import MosfetSimulation
from .simulators.spice_simulators.utils import sweep_axis
from .table_cleanup import LookupTableCleaner

# >>>

class LookupTableGenerator:
    def __init__(
        self,
        *,
        simulator,
        model_sweeps,
        n_process=1,
        description="gmid lookup table",
    ):
        self.model_sweeps = model_sweeps
        self.simulator = simulator
        self.description = description
        self.n_process = n_process

    def op_simulation(self):
        simulation = MosfetSimulation(
            self.simulator,
            self.model_sweeps,
            self.n_process
        )
        simulation.op_simulation()

    def build(self, filepath):
        # Create MosfetSimulation instance.
        simulation = MosfetSimulation(
            self.simulator,
            self.model_sweeps,
            self.n_process
        )

        # Run all simulation jobs with progress updates.
        simulation.simulate()
        lookup_table = simulation.lookup_table
        parameters_to_save = [key.lower() for key in self.simulator.parameters_to_save]
        device_parameters = {key.lower(): value for key, value in self.simulator.device_parameters.items()}

        # Cleanup table to remove entries for parameters that were not found
        # or parameters that have a constant value throughout.
        cleaner = LookupTableCleaner(lookup_table, parameters_to_save)
        cleaner.clean_lookup_table()

        # Store general and grid information.
        lookup_table["description"] = self.description
        lookup_table["simulator"] = self.simulator.__class__.__name__
        lookup_table["parameter_names"] = parameters_to_save
        lookup_table["device_parameters"] = device_parameters
        for transistor_name, sweep in self.model_sweeps.items():
            per_transistor_params = dict(device_parameters)
            per_transistor_params.update(cleaner.scalar_params.get(transistor_name, {}))
            lookup_table[transistor_name]["vgs"] = sweep_axis(sweep.vgs)
            lookup_table[transistor_name]["vds"] = sweep_axis(sweep.vds)
            lookup_table[transistor_name]["vbs"] = sweep_axis(sweep.vbs)
            lookup_table[transistor_name]["length"] = np.array(sweep.length)
            lookup_table[transistor_name]["model_name"] = transistor_name
            lookup_table[transistor_name]["parameter_names"] = parameters_to_save
            lookup_table[transistor_name]["device_parameters"] = per_transistor_params

        directory = os.path.dirname(filepath)
        if directory and not os.path.exists(directory):
            os.makedirs(directory)

        np.savez_compressed(f"{filepath}.npz", lookup_table=np.array(lookup_table, dtype=object))
        print("Done")
