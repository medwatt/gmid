# imports <<<
import os
import pickle
import tempfile
import subprocess
import numpy as np
from .base_simulator import BaseSimulator
from .parsers.hspice import import_export
# >>>

class HspiceSimulator(BaseSimulator):
    def __init__(
        self,
        simulator_path,
        temperature,
        model_paths,
        parameters_to_save,
        mos_spice_symbols,
    ):
        super().__init__(simulator_path, temperature, model_paths, parameters_to_save)
        self.mos_spice_symbols = mos_spice_symbols
        self.make_temp_files()

    def make_temp_files(self):
        self.tmp_dir = tempfile.mkdtemp()
        self.input_file_path = os.path.join(self.tmp_dir, "input.txt")
        self.log_file_path = os.path.join(self.tmp_dir, "input.lis")
        self.output_file_path = os.path.join(self.tmp_dir, "input.sw0")
        self.decoded_output_path = os.path.join(self.tmp_dir, "input_sw0.pickle")

    def setup_op_simulation(self, vgs, vds):
        return [
            f".TEMP = {self.temperature}",
            ".option POST=2",
            ".op",
            ".end",
        ]

    def setup_dc_simulation(self, vgs, vds):
        vgs_start, vgs_stop, vgs_step = vgs
        vds_start, vds_stop, vds_step = vds
        analysis_string = f".dc VGS {vgs_start} {vgs_stop} {vgs_step} VDS {vds_start} {vds_stop} {vds_step}"

        symbol = self.mos_spice_symbols[1]
        self.parameter_table = {
            # parameter name : [name recognized by simulator, name used in the output file],
            "id": [
                ".probe DC m_id = par('abs(i(vds))')",
                "m_id"
            ],
            "vth": [
                f".probe DC m_vth = par('vth({symbol})')",
                "m_vth"
            ],
            "vdsat": [
                f".probe DC m_vdsat = par('vdsat({symbol})')",
                "m_vdsat",
            ],
            "gm": [
                f".probe DC m_gm = par('gmo({symbol})')",
                "m_gm"
            ],
            "gmbs": [
                f".probe DC m_gmb = par('gmbso({symbol})')",
                "m_gmb",
            ],
            "gds": [
                f".probe DC m_gds = par('gdso({symbol})')",
                "m_gds",
            ],
            "cgg": [
                f".probe DC m_cgg = par('cggbo({symbol})')",
                "m_cgg",
            ],
            "cgs": [
                f".probe DC m_cgs = par('-cgsbo({symbol})')",
                "m_cgs",
            ],
            "cgd": [
                f".probe DC m_cgd = par('-cgdbo({symbol})')",
                "m_cgd",
            ],
            "cgb": [
                f".probe DC m_cgb = par('cggbo({symbol})-(-cgsbo({symbol}))-(-cgdbo({symbol}))')",
                "m_cgb",
            ],
            "cdd": [
                f".probe DC m_cdd = par('cddbo({symbol})')",
                "m_cdd",
            ],
            "css": [
                f".probe DC m_css = par('-cgsbo({symbol})-cbsbo({symbol})')",
                "m_css",
            ],
        }

        self.parameter_table = {
            k: v
            for k, v in self.parameter_table.items()
            if k in self.parameters_to_save
        }

        return [
            f".TEMP = {self.temperature}",
            ".options probe dccap brief accurate",
            ".option POST=2",
            "\n".join([val[0] for val in self.parameter_table.values()]),
            analysis_string,
            ".end",
        ]

    def run_simulation(self, netlist, verbose=False):
        with open(self.input_file_path, "w") as f:
            f.write("\n".join(netlist))

        if verbose:
            cmd = f"{self.simulator_path} {self.input_file_path}"
            with subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, bufsize=1, universal_newlines=True) as p:
                for line in p.stdout:
                    print(line, end='')

        else:
            cmd = f"{self.simulator_path} -i {self.input_file_path} -o {self.tmp_dir}"
            result = subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

            if result.returncode != 0:
                with open(self.log_file_path, "r") as log_file:
                    log_contents = log_file.read()
                self.remove_temp_files()
                raise Exception(f"Simulation error: {log_contents}")

    def parse_output(self):
        import_export(self.output_file_path, "pickle")
        with open(self.decoded_output_path, "rb") as file:
            loaded_data = pickle.load(file)
        return loaded_data

    def save_parameters(self, analysis, transistor_type, length, vbs, lookup_table, n_vgs, n_vds):
        for p in self.parameters_to_save:
            col_name = self.parameter_table[p][1]
            if col_name in analysis.keys():
                res = np.array(analysis[col_name]).T
                lookup_table[transistor_type][p][length][vbs] = res
