# imports <<<
import os
import tempfile

import numpy as np

from .base_simulator import BaseSimulator
from .parsers.spectre import parse_file
# >>>

class SpectreSimulator(BaseSimulator):
    def __init__(
        self,
        *,
        temperature=27,
        raw_spice=None,
        lib_mappings=None,
        include_paths=None,
        simulator_path="spectre",
        mos_spice_symbols=("M1", "M1"),
        device_parameters=None,
        parameters_to_save=None,
    ):
        super().__init__(
            raw_spice=raw_spice,
            temperature=temperature,
            lib_mappings=lib_mappings,
            include_paths=include_paths,
            simulator_path=simulator_path,
            mos_spice_symbols=mos_spice_symbols,
            device_parameters=device_parameters if device_parameters is not None else {"w": 10e-6},
            parameters_to_save=parameters_to_save if parameters_to_save is not None
            else ["id", "vth", "vdsat", "gm", "gmbs", "gds", "cgg", "cgs", "cgb", "cgd", "cdd"],
        )

    def make_temp_files(self):
        self.tmp_dir = tempfile.mkdtemp()
        self.input_file_path = os.path.join(self.tmp_dir, "input.scs")
        self.log_file_path = os.path.join(self.tmp_dir, "spectre.log")
        self.output_file_path = os.path.join(self.tmp_dir, "output.raw")

    def build_simulation_command(self, verbose):
        if verbose:
            return [self.simulator_path, "-format", "nutbin", "-raw", self.output_file_path, self.input_file_path]
        return [self.simulator_path, "-format", "nutbin", "-raw", self.output_file_path, "=log", self.log_file_path, self.input_file_path]

    def setup_op_simulation(self, sweep):
        symbol = self.mos_spice_symbols[1]
        return [
            f"o_temp options temp={self.temperature} tnom={self.temperature}",
            f"save {symbol}:oppoint",
            "dcop dc oppoint=screen",
        ]

    def setup_dc_simulation(self, sweep):
        symbol = self.mos_spice_symbols[1]
        id_terminal = "p" if sweep.mos_type == "pmos" else "n"
        id_probe = f"VDS:{id_terminal}"

        self.parameter_table = {
            "id":    id_probe,
            "vth":   f"{symbol}:vth",
            "vdsat": (f"{symbol}:vdsat", f"{symbol}:vdss"),  # BSIM: vdsat; PSP: vdss
            "gm":    f"{symbol}:gm",
            "gmbs":  (f"{symbol}:gmbs", f"{symbol}:gmb"),  # BSIM: gmbs; PSP: gmb
            "gds":   f"{symbol}:gds",
            "cgg":   f"{symbol}:cgg",
            "cgs":   f"{symbol}:cgs",
            "cgd":   f"{symbol}:cgd",
            "cgb":   f"{symbol}:cgb",
            "cdd":   f"{symbol}:cdd",
            "css":   f"{symbol}:css",
        }

        self.parameter_table = self.select_parameters(self.parameter_table)

        vgs_start, vgs_stop, vgs_step = sweep.vgs
        vds_start, vds_stop, vds_step = sweep.vds

        return [
            f"o_temp options temp={self.temperature} tnom={self.temperature}",
            f"save {symbol}:oppoint {id_probe}",
            f"vgs_swp sweep param=dc dev=VGS start={vgs_start} stop={vgs_stop} step={vgs_step} {{",
            f"    vds_dc dc param=dc dev=VDS start={vds_start} stop={vds_stop} step={vds_step}",
            "}",
        ]

    def parse_output(self):
        arrs, _ = parse_file(self.output_file_path)
        return arrs

    def extract_parameters(self, analysis, n_vgs, n_vds):
        results = {}
        for p, names in self.parameter_table.items():
            # a model may name an output differently: the first name present wins
            for col_name in (names,) if isinstance(names, str) else names:
                try:
                    if len(analysis) == 1:
                        data = np.array(analysis[0][col_name]).reshape(n_vgs, n_vds)
                    else:
                        data = np.array([arr[col_name] for arr in analysis])
                    results[p] = data.astype(np.float32)
                    break
                except (KeyError, ValueError):
                    pass
        return results

