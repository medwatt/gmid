# ./lookup_table_generator.py
# imports <<<
import os
import numpy as np
from .simulators.ngspice_simulator import NgspiceSimulator
from .simulators.hspice_simulator import HspiceSimulator
from .mosfet_netlist_generator import MosfetNetlistGenerator
from .mosfet_simulation import MosfetSimulation
# >>>

class LookupTableGenerator:
    def __init__(
        self,
        *,
        model_sweeps,                       # model_sweeps now = { "NMOS_VTH": TransistorSweep(...), ...}
        include_paths=None,                 # include_paths = [ "./models/NMOS_VTH.lib", ...]
        lib_mappings=None,                  # lib_mappings = [("./models/design_wrapper.lib", "tt_pre") ...],
        raw_spice=None,                     # raw_spice = ["", ...]
        mos_spice_symbols=("m1", "m1"),     # mos_spice_symbols = ("mosfet/subcircuit name", "hierarchical_name_of_transistor")
        width=10e-6,
        simulator="ngspice",
        simulator_path=None,
        temperature=25,
        parameters_to_save=[],
        description="gmid lookup table",
    ):
        self.model_sweeps = model_sweeps
        self.include_paths = include_paths
        self.lib_mappings = lib_mappings
        self.raw_spice = raw_spice
        self.mos_spice_symbols = mos_spice_symbols
        self.width = width
        self.simulator = simulator
        self.simulator_path = simulator_path or simulator
        self.temperature = temperature
        self.parameters_to_save = parameters_to_save or ["id", "vth", "vdsat", "gm", "gmbs", "gds", "cgg", "cgs", "cgb", "cgd", "cdd"]
        self.description = description
        self.lookup_table = {}

    def choose_simulator(self):
        if self.simulator.lower() == "ngspice":
            return NgspiceSimulator(
                self.simulator_path,
                self.temperature,
                self.include_paths,
                self.parameters_to_save,
                self.mos_spice_symbols,
            )
        elif self.simulator.lower() == "hspice":
            return HspiceSimulator(
                self.simulator_path,
                self.temperature,
                self.include_paths,
                self.parameters_to_save,
                self.mos_spice_symbols,
            )
        else:
            raise ValueError("Unsupported simulator type.")

    def prepare_simulation(self):
        simulator = self.choose_simulator()
        netlist_gen = MosfetNetlistGenerator(
            self.model_sweeps,
            self.width,
            self.mos_spice_symbols,
            self.include_paths,
            self.lib_mappings,
            self.raw_spice,
        )
        return MosfetSimulation(simulator, netlist_gen, self.model_sweeps)

    def op_simulation(self):
        simulator = self.choose_simulator()
        simulation = self.prepare_simulation()
        simulation.op_simulation()

    def build(self, filepath):
        def range_args(tup):
            start, stop, step = tup
            return (start, stop + step, step)

        simulator = self.choose_simulator()
        simulation = self.prepare_simulation()
        simulation.simulate()
        self.lookup_table = simulation.lookup_table

        # General table information
        self.lookup_table["width"] = self.width
        self.lookup_table["description"] = self.description
        self.lookup_table["simulator"] = self.simulator
        self.lookup_table["parameter_names"] = self.parameters_to_save

        # Store grid and meta-data for each transistor model
        for transistor_name, sweep in self.model_sweeps.items():
            self.lookup_table[transistor_name]["vgs"] = np.arange(*range_args(sweep.vgs))
            self.lookup_table[transistor_name]["vds"] = np.arange(*range_args(sweep.vds))
            self.lookup_table[transistor_name]["vbs"] = np.arange(*range_args(sweep.vbs))
            self.lookup_table[transistor_name]["length"] = np.array(sweep.length)
            self.lookup_table[transistor_name]["width"] = self.width
            self.lookup_table[transistor_name]["model_name"] = transistor_name
            self.lookup_table[transistor_name]["parameter_names"] = self.parameters_to_save

        # Ensure the directory exists
        directory = os.path.dirname(filepath)
        if directory and not os.path.exists(directory):
            os.makedirs(directory)

        # Save the file
        np.savez_compressed(f"{filepath}.npz", lookup_table=np.array(self.lookup_table, dtype=object))

        # Clean up
        simulator.remove_temp_files()
        print("Done")
