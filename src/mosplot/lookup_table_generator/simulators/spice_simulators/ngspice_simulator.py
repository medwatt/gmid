# imports <<<
import os
import tempfile

import numpy as np

from .base_simulator import BaseSimulator
from .parsers.ngspice import NgspiceRawFileReader
# >>>

class NgspiceSimulator(BaseSimulator):
    def __init__(
        self,
        *,
        temperature=27,
        raw_spice=None,
        lib_mappings=None,
        include_paths=None,
        simulator_path="ngspice",
        mos_spice_symbols=("m1", "m1"),
        parameters_to_save=["id", "vth", "vdsat", "gm", "gmbs", "gds", "cgg", "cgs", "cgb", "cgd", "cdd"],
    ):
        super().__init__(
                raw_spice=raw_spice,
                temperature=temperature,
                lib_mappings=lib_mappings,
                include_paths=include_paths,
                simulator_path=simulator_path,
                mos_spice_symbols=mos_spice_symbols,
                parameters_to_save=parameters_to_save,
        )

    def make_temp_files(self):
        self.tmp_dir = tempfile.mkdtemp()
        self.input_file_path = os.path.join(self.tmp_dir, "input.txt")
        self.log_file_path = os.path.join(self.tmp_dir, "log.txt")
        self.output_file_path = os.path.join(self.tmp_dir, "output.txt")

    def build_simulation_command(self, verbose):
        if verbose:
            return f"{self.simulator_path} -b {self.input_file_path}"
        return f"{self.simulator_path} -b -o {self.log_file_path} {self.input_file_path}"

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
        vgs_start, vgs_stop, vgs_step = vgs
        vds_start, vds_stop, vds_step = vds
        analysis_string = f"dc VDS {vds_start} {vds_stop} {vds_step} VGS {vgs_start} {vgs_stop} {vgs_step}"
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

    def parse_output(self):
        reader = NgspiceRawFileReader()
        analysis, _ = reader.read_file(self.output_file_path)
        return analysis

    def extract_parameters(self, analysis, n_vgs, n_vds):
        results = {}
        column_names = analysis[0].dtype.names
        data = analysis[0]
        for p in self.parameter_table.keys():
            col_name = self.parameter_table[p][1]
            if col_name in column_names:
                res = np.array(data[col_name]).reshape(n_vgs, n_vds)
                results[p] = res
        return results

    def save_parameters(self, analysis, transistor_type, length, vbs, lookup_table, n_vgs, n_vds):
        column_names = analysis[0].dtype.names
        data = analysis[0]
        for p in self.parameter_table.keys():
            col_name = self.parameter_table[p][1]
            if col_name in column_names:
                res = np.array(data[col_name]).reshape(n_vgs, n_vds)
                lookup_table[transistor_type][p][length][vbs] = res
