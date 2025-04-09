# imports <<<
import os
import tempfile
import subprocess
import numpy as np
from .base_simulator import BaseSimulator
from .parsers.ngspice import NgspiceRawFileReader
from .mosfet_simulation import MosfetSimulation
from .spice_mosfet_netlist_generator import SpiceMosfetNetlistGenerator
# >>>

class NgspiceSimulator(BaseSimulator):
    def __init__(
        self,
        *,
        temperature=27,
        raw_spice=None,
        lib_mappings=None,
        include_paths=None,
        simulator_path=None,
        mos_spice_symbols=("m1", "m1"),
        parameters_to_save=["id", "vth", "vdsat", "gm", "gmbs", "gds", "cgg", "cgs", "cgb", "cgd", "cdd"],
    ):
        super().__init__(simulator_path, temperature, include_paths, lib_mappings, parameters_to_save)
        self.mos_spice_symbols = mos_spice_symbols
        self.raw_spice = raw_spice
        self.make_temp_files()

    def make_temp_files(self):
        self.tmp_dir = tempfile.mkdtemp()
        self.input_file_path = os.path.join(self.tmp_dir, "input.txt")
        self.log_file_path = os.path.join(self.tmp_dir, "log.txt")
        self.output_file_path = os.path.join(self.tmp_dir, "output.txt")

    def setup_op_simulation(self, vgs, vds):
        return [
            f".options TEMP = {self.temperature}",
            f".options TNOM = {self.temperature}",
            ".control",
            "op",
            "show all",
            ".endc",
            ".end",
        ]

    def setup_dc_simulation(self, vgs, vds):
        vgs_start, vgs_stop, vgs_step = vgs
        vds_start, vds_stop, vds_step = vds
        analysis_string = f"dc VDS {vds_start} {vds_stop} {vds_step} VGS {vgs_start} {vgs_stop} {vgs_step}"

        symbol = self.mos_spice_symbols[1]
        self.parameter_table = {
            "id":    ["save i(vds)",            "i(i_vds)"],
            "vth":   [f"save @{symbol}[vth]",   f"v(@{symbol}[vth])"],
            "vdsat": [f"save @{symbol}[vdsat]", f"v(@{symbol}[vdsat])"],
            "gm":    [f"save @{symbol}[gm]",    f"@{symbol}[gm]"],
            "gmbs":  [f"save @{symbol}[gmbs]",  f"@{symbol}[gmbs]"],
            "gds":   [f"save @{symbol}[gds]",   f"@{symbol}[gds]"],
            "cgg":   [f"save @{symbol}[cgg]",   f"@{symbol}[cgg]"],
            "cgs":   [f"save @{symbol}[cgs]",   f"@{symbol}[cgs]"],
            "cbg":   [f"save @{symbol}[cbg]",   f"@{symbol}[cbg]"],
            "cgd":   [f"save @{symbol}[cgd]",   f"@{symbol}[cgd]"],
            "cdd":   [f"save @{symbol}[cdd]",   f"@{symbol}[cdd]"],
        }

        self.parameter_table = { k: v for k, v in self.parameter_table.items() if k in self.parameters_to_save }

        return [
            f".options TEMP = {self.temperature}",
            f".options TNOM = {self.temperature}",
            ".control",
            "\n".join([val[0] for val in self.parameter_table.values()]),
            analysis_string,
            "let i_vds = abs(i(vds))",
            f"write {self.output_file_path} all",
            ".endc",
            ".end",
        ]

    def run_simulation(self, netlist, verbose=False):
        with open(self.input_file_path, "w") as f:
            f.write("\n".join(netlist))

        if verbose:
            cmd = f"{self.simulator_path} -b {self.input_file_path}"
            with subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                  bufsize=1, universal_newlines=True) as p:
                for line in p.stdout:
                    print(line, end='')
        else:
            cmd = f"{self.simulator_path} -b -o {self.log_file_path} {self.input_file_path}"
            result = subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            if result.returncode != 0:
                with open(self.log_file_path, "r") as log_file:
                    log_contents = log_file.read()
                self.remove_temp_files()
                raise Exception(f"Simulation error: {log_contents}")

    def parse_output(self):
        reader = NgspiceRawFileReader()
        analysis, _ = reader.read_file(self.output_file_path)
        return analysis

    def save_parameters(self, analysis, transistor_type, length, vbs, lookup_table, n_vgs, n_vds):
        column_names = analysis[0].dtype.names
        data = analysis[0]
        for p in self.parameter_table.keys():
            col_name = self.parameter_table[p][1]
            if col_name in column_names:
                res = np.array(data[col_name])
                lookup_table[transistor_type][p][length][vbs] = res.reshape(n_vgs, n_vds)

    def prepare_simulation(self, model_sweeps, width):
        netlist_gen = SpiceMosfetNetlistGenerator(
            model_sweeps,
            width,
            self.mos_spice_symbols,
            self.include_paths,
            self.lib_mappings,
            self.raw_spice,
        )
        return MosfetSimulation(self, netlist_gen, model_sweeps)
