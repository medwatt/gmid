# imports <<<
import os
import tempfile

import numpy as np

from .base_simulator import BaseSimulator
from .parsers.ngspice import parse_file
# >>>

class NgspiceSimulator(BaseSimulator):
    def __init__(
        self,
        *,
        temperature=27,
        raw_spice=None,
        lib_mappings=None,
        include_paths=None,
        osdi_paths=None,
        simulator_path="ngspice",
        mos_spice_symbols=("m1", "m1"),
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
                else ["weff", "id", "vth", "vdsat", "vdssat", "gm", "gmbs", "gds", "cgg", "cgs", "cbg", "cgd", "cdd"],
        )
        self.osdi_paths = osdi_paths
        self._init_config["osdi_paths"] = osdi_paths

    def make_temp_files(self):
        self.tmp_dir = tempfile.mkdtemp()
        self.input_file_path = os.path.join(self.tmp_dir, "input.txt")
        self.log_file_path = os.path.join(self.tmp_dir, "log.txt")
        self.output_file_path = os.path.join(self.tmp_dir, "output.txt")

    def build_simulation_command(self, verbose):
        if verbose:
            return [self.simulator_path, "-b", self.input_file_path]
        return [self.simulator_path, "-b", "-o", self.log_file_path, self.input_file_path]

    def setup_op_simulation(self, sweep):
        osdi = None
        if self.osdi_paths:
            osdi = "\n".join([f"pre_osdi {p}" for p in self.osdi_paths])
        return [
            f".options TEMP = {self.temperature}",
            f".options TNOM = {self.temperature}",
            ".control",
            osdi,
            "op",
            "show all",
            ".endc",
            ".end",
        ]

    def setup_dc_simulation(self, sweep):
        symbol = self.mos_spice_symbols[1]
        v_sign = "-" if sweep.mos_type == "pmos" else ""
        self.parameter_table = {
            "id":     ["save i(vds)",             "i(i_vds)"],
            "weff":   [f"save @{symbol}[weff]",   f"v(@{symbol}[weff])"],
            "vth":    [f"save @{symbol}[vth]",    "v(m_vth)"],
            "vdsat":  [f"save @{symbol}[vdsat] @{symbol}[vdss]", ("v(m_vdsat)", "v(m_vdss)")],  # BSIM: vdsat; PSP: vdss
            "vdssat": [f"save @{symbol}[vdssat]", "v(m_vdssat)"],
            "gm":     [f"save @{symbol}[gm]",     f"@{symbol}[gm]"],
            "gmbs":   [f"save @{symbol}[gmbs] @{symbol}[gmb]", (f"@{symbol}[gmbs]", f"@{symbol}[gmb]")],  # BSIM: gmbs; PSP: gmb
            "gds":    [f"save @{symbol}[gds]",    f"@{symbol}[gds]"],
            "cgg":    [f"save @{symbol}[cgg]",    f"@{symbol}[cgg]"],
            "cgs":    [f"save @{symbol}[cgs]",    f"@{symbol}[cgs]"],
            "cbg":    [f"save @{symbol}[cbg]",    f"@{symbol}[cbg]"],
            "cgd":    [f"save @{symbol}[cgd]",    f"@{symbol}[cgd]"],
            "cdd":    [f"save @{symbol}[cdd]",    f"@{symbol}[cdd]"],
        }
        self.parameter_table = self.select_parameters(self.parameter_table)
        vgs_start, vgs_stop, vgs_step = sweep.vgs
        vds_start, vds_stop, vds_step = sweep.vds
        analysis_string = f"dc VDS {vds_start} {vds_stop} {vds_step} VGS {vgs_start} {vgs_stop} {vgs_step}"
        osdi = None
        if self.osdi_paths:
            osdi = "\n".join([f"pre_osdi {p}" for p in self.osdi_paths])

        # Derived vectors are only emitted for requested parameters. A `let` on a
        # device parameter the model does not expose (vdssat on a non-BSIM4 model,
        # vdsat on PSP) prints an error; ngspice skips it and the run continues.
        derived_vectors = {
            "id": "let i_vds = abs(i(vds))",
            "vth": f"let m_vth = {v_sign}abs(@{symbol}[vth])",
            "vdsat": f"let m_vdsat = {v_sign}abs(@{symbol}[vdsat])\nlet m_vdss = {v_sign}abs(@{symbol}[vdss])",
            "vdssat": f"let m_vdssat = {v_sign}abs(@{symbol}[vdssat])",
        }
        let_lines = [line for p, line in derived_vectors.items() if p in self.parameter_table]

        return [
            f".options TEMP = {self.temperature}",
            f".options TNOM = {self.temperature}",
            ".control",
            osdi,
            "\n".join([val[0] for val in self.parameter_table.values()]),
            analysis_string,
            *let_lines,
            f"write {self.output_file_path} all",
            ".endc",
            ".end",
        ]

    def parse_output(self):
        analysis, _ = parse_file(self.output_file_path)
        return analysis

    def extract_parameters(self, analysis, n_vgs, n_vds):
        results = {}
        column_names = analysis[0].dtype.names
        data = analysis[0]
        for p, (_, names) in self.parameter_table.items():
            # a model may name an output differently: the first name present wins
            for col_name in (names,) if isinstance(names, str) else names:
                if col_name in column_names:
                    results[p] = np.array(data[col_name]).reshape(n_vgs, n_vds)
                    break
        return results
