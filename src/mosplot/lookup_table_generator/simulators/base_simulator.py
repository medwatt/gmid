# import <<<
from abc import ABC, abstractmethod
import os
import shutil
# >>>

class BaseSimulator(ABC):
    def __init__(self, simulator_path, temperature, model_paths, parameters_to_save):
        self.simulator_path = simulator_path
        self.temperature = temperature
        self.model_paths = model_paths
        self.parameters_to_save = parameters_to_save
        self.tmp_dir = None
        self.input_file_path = None
        self.log_file_path = None
        self.output_file_path = None
        self.parameter_table = {}

    def validate_paths(self):
        paths = self.model_paths.copy()
        if self.simulator_path:
            paths.append(self.simulator_path)
        for path in paths:
            if not os.path.exists(path):
                raise FileNotFoundError(f"Path {path} not found.")

    def setup(self):
        self.validate_paths()
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
    def save_parameters(self, analysis, transistor_type, length, vsb, lookup_table, n_vgs, n_vds):
        pass
