import os
import numpy as np
from .simulators.spice_simulators.mosfet_simulation import MosfetSimulation

class LookupTableGenerator:
    def __init__(
        self,
        *,
        simulator,
        model_sweeps,
        width=10e-6,
        n_process=1,
        description="gmid lookup table",
    ):
        self.model_sweeps = model_sweeps
        self.simulator = simulator
        self.width = width
        self.description = description
        self.n_process = n_process
        self.lookup_table = {}

    def op_simulation(self):
        simulation = MosfetSimulation(
            self.simulator,
            self.model_sweeps,
            self.width,
            self.n_process
        )

        simulation.print_netlist()
        simulation.op_simulation()

    def build(self, filepath):
        def range_args(tup):
            start, stop, step = tup
            return (start, stop + step, step)

        # Create MosfetSimulation instance.
        simulation = MosfetSimulation(
            self.simulator,
            self.model_sweeps,
            self.width,
            self.n_process
        )

        # Print a sample netlist before simulation starts.
        simulation.print_netlist()

        # Run all simulation jobs with progress updates.
        simulation.simulate()
        self.lookup_table = simulation.lookup_table

        # Store general and grid information.
        self.lookup_table["width"] = self.width
        self.lookup_table["description"] = self.description
        self.lookup_table["simulator"] = self.simulator.__class__.__name__
        self.lookup_table["parameter_names"] = self.simulator.parameters_to_save
        for transistor_name, sweep in self.model_sweeps.items():
            self.lookup_table[transistor_name]["vgs"] = np.arange(*range_args(sweep.vgs))
            self.lookup_table[transistor_name]["vds"] = np.arange(*range_args(sweep.vds))
            self.lookup_table[transistor_name]["vbs"] = np.arange(*range_args(sweep.vbs))
            self.lookup_table[transistor_name]["length"] = np.array(sweep.length)
            self.lookup_table[transistor_name]["width"] = self.width
            self.lookup_table[transistor_name]["model_name"] = transistor_name
            self.lookup_table[transistor_name]["parameter_names"] = self.simulator.parameters_to_save

        directory = os.path.dirname(filepath)
        if directory and not os.path.exists(directory):
            os.makedirs(directory)

        np.savez_compressed(f"{filepath}.npz", lookup_table=np.array(self.lookup_table, dtype=object))
        self.simulator.remove_temp_files()
        print("Done")
