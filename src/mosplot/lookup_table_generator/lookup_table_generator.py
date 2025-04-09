# imports <<<
import os
import numpy as np
# >>>

class LookupTableGenerator:
    def __init__(
        self,
        *,
        model_sweeps,
        simulator,
        width=10e-6,
        description="gmid lookup table",
    ):
        self.model_sweeps = model_sweeps
        self.simulator = simulator
        self.width = width
        self.description = description
        self.lookup_table = {}

    def op_simulation(self):
        simulation = self.simulator.prepare_simulation(self.model_sweeps, self.width)
        simulation.op_simulation()

    def build(self, filepath):
        def range_args(tup):
            start, stop, step = tup
            return (start, stop + step, step)

        simulation = self.simulator.prepare_simulation(self.model_sweeps, self.width)
        simulation.simulate()
        self.lookup_table = simulation.lookup_table

        # General table information.
        self.lookup_table["width"] = self.width
        self.lookup_table["description"] = self.description
        self.lookup_table["simulator"] = self.simulator.__class__.__name__
        self.lookup_table["parameter_names"] = self.simulator.parameters_to_save

        # Store grid and meta-data for each transistor model.
        for transistor_name, sweep in self.model_sweeps.items():
            self.lookup_table[transistor_name]["vgs"] = np.arange(*range_args(sweep.vgs))
            self.lookup_table[transistor_name]["vds"] = np.arange(*range_args(sweep.vds))
            self.lookup_table[transistor_name]["vbs"] = np.arange(*range_args(sweep.vbs))
            self.lookup_table[transistor_name]["length"] = np.array(sweep.length)
            self.lookup_table[transistor_name]["width"] = self.width
            self.lookup_table[transistor_name]["model_name"] = transistor_name
            self.lookup_table[transistor_name]["parameter_names"] = self.simulator.parameters_to_save

        # Ensure the directory exists.
        directory = os.path.dirname(filepath)
        if directory and not os.path.exists(directory):
            os.makedirs(directory)

        # Save the file.
        np.savez_compressed(f"{filepath}.npz", lookup_table=np.array(self.lookup_table, dtype=object))

        # Clean up temporary simulator files.
        self.simulator.remove_temp_files()
        print("Done")
