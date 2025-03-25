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
        lengths,                            # lengths = [45e-9, 100e-9, 200e-9, ...],
        model_names,                        # model_names = {"NMOS_VTH": "nmos", "PMOS_VTH": "pmos", ...}
        include_paths=None,                 # include_paths = [ "./models/NMOS_VTH.lib", ...]
        lib_path_names=None,                # lib_path_names = [("./models/design_wrapper.lib", tt_pre") ...],
        raw_spice=None,                     # raw_spice = ["", ...]
        mos_spice_symbols=("m1", "m1"),     # mos_spice_symbols = ("mosfet/subcircuit name", "hierarchical_name_of_transistor")
        width=10e-6,
        vgs=(0, 1, 0.01),
        vds=(0, 1, 0.01),
        vsb=(0, 1, 0.1),
        simulator="ngspice",
        simulator_path=None,
        temperature=27,
        parameters_to_save=[],
        description="gmid lookup table",
    ):
        self.lengths = np.array(lengths)
        self.model_names = model_names
        self.include_paths = include_paths
        self.lib_path_names = lib_path_names
        self.raw_spice = raw_spice
        self.mos_spice_symbols = mos_spice_symbols
        self.width = width
        self.vgs = np.array(vgs)
        self.vds = np.array(vds)
        self.vsb = np.array(vsb)
        self.simulator_type = simulator
        self.simulator_path = simulator_path or simulator
        self.temperature = temperature
        self.parameters_to_save = parameters_to_save
        self.description = description
        self.lookup_table = {}

    def _range_args(self, tup):
        start, stop, step = tup
        return (start, stop + step, step)

    def choose_simulator(self):
        if self.simulator_type.lower() == "ngspice":
            if not self.parameters_to_save:
                self.parameters_to_save = ["id", "vth", "vdsat", "gm", "gmbs", "gds", "cgg", "cgs", "cbg", "cgd", "cdd"]
            return NgspiceSimulator(
                self.simulator_path,
                self.temperature,
                self.include_paths,
                self.parameters_to_save,
                self.mos_spice_symbols,
            )
        elif self.simulator_type.lower() == "hspice":
            if not self.parameters_to_save:
                self.parameters_to_save = ["id", "vth", "vdsat", "gm", "gmbs", "gds", "cgg", "cgs", "cgb", "cgd", "cdd", "css"]
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
            self.model_names,
            self.width,
            self.mos_spice_symbols,
            self.include_paths,
            self.lib_path_names,
            self.raw_spice,
        )

        n_vgs = int(round((self.vgs[1] - self.vgs[0]) / self.vgs[2])) + 1
        n_vds = int(round((self.vds[1] - self.vds[0]) / self.vds[2])) + 1
        n_vsb = int(round((self.vsb[1] - self.vsb[0]) / self.vsb[2])) + 1

        return MosfetSimulation(
            simulator,
            netlist_gen,
            self.vgs,
            self.vds,
            self.vsb,
            self.lengths,
            n_vgs,
            n_vds,
            n_vsb,
            self.model_names,
        )

    def op_simulation(self):
        simulator = self.choose_simulator()
        simulation = self.prepare_simulation()
        simulation.op_simulation()

    def build(self, filepath):
        simulator = self.choose_simulator()
        simulation = self.prepare_simulation()
        simulation.simulate()
        self.lookup_table = simulation.lookup_table

        # General table information
        self.lookup_table["description"] = self.description
        self.lookup_table["parameter_names"] = self.parameters_to_save
        self.lookup_table["width"] = self.width
        self.lookup_table["lengths"] = self.lengths

        # Store grid and meta-data for each transistor model
        for transistor_name, transistor_type in self.model_names.items():
            r = -1 if transistor_type == "pmos" else 1
            self.lookup_table[transistor_name]["vgs"] = np.arange(*self._range_args(self.vgs*r))
            self.lookup_table[transistor_name]["vds"] = np.arange(*self._range_args(self.vds*r))
            self.lookup_table[transistor_name]["vsb"] = np.arange(*self._range_args(self.vsb*r))
            self.lookup_table[transistor_name]["width"] = self.width
            self.lookup_table[transistor_name]["lengths"] = self.lengths
            self.lookup_table[transistor_name]["model_name"] = transistor_name
            self.lookup_table[transistor_name]["model_type"] = transistor_type
            self.lookup_table[transistor_name]["parameter_names"] = self.parameters_to_save

        # Ensure the directory exists
        directory = os.path.dirname(filepath)
        if directory and not os.path.exists(directory):
            os.makedirs(directory)

        # Save the file
        np.save(f"{filepath}.npy", self.lookup_table, allow_pickle=True)

        simulator.remove_temp_files()
        print("Done")
