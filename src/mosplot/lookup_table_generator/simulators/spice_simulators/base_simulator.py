# imports <<<
from abc import ABC, abstractmethod
import os
import shutil
# >>>

class BaseSimulator(ABC):
    def __init__(self, simulator_path, temperature, include_paths, lib_mappings, parameters_to_save):
        self.simulator_path = simulator_path
        self.temperature = temperature
        self.include_paths = include_paths
        self.lib_mappings = lib_mappings
        self.parameters_to_save = parameters_to_save
        self.tmp_dir = None
        self.validate_paths()

    def validate_paths(self):
        paths = []

        # Check each include path.
        if self.include_paths:
            for p in self.include_paths:
                paths.append(p)

        if self.lib_mappings:
            # Check first entry of every lib mapping tuple.
            for item in self.lib_mappings:
                paths.append(item[0])

        for path in paths:
            if not os.path.exists(path):
                raise FileNotFoundError(f"Path {path} not found.")

        # Make sure the simulator binary exists.
        if not shutil.which(self.simulator_path):
            raise ValueError(f"Binary '{self.simulator_path}' not accessible.")

    @abstractmethod
    def make_temp_files(self):
        pass

    def remove_temp_files(self):
        if self.tmp_dir and os.path.isdir(self.tmp_dir):
            shutil.rmtree(self.tmp_dir)

    @abstractmethod
    def setup_op_simulation(self, vgs, vds):
        pass

    @abstractmethod
    def setup_dc_simulation(self, vgs, vds):
        pass

    @abstractmethod
    def run_simulation(self, netlist, verbose):
        pass

    @abstractmethod
    def parse_output(self):
        pass

    @abstractmethod
    def save_parameters(self, analysis, transistor_type, length, vbs, lookup_table, n_vgs, n_vds):
        pass
