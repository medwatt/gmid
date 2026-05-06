class SpectreMosfetNetlistGenerator:
    def __init__(self, model_sweeps, device_parameters, mos_spice_symbols, include_paths, lib_path_and_names, raw_spice):
        self.model_sweeps = model_sweeps
        self.raw_spice = raw_spice
        self.device_parameters = device_parameters
        self.mos_spice_symbols = mos_spice_symbols
        self.include_paths = include_paths
        self.lib_path_and_names = lib_path_and_names

    def get_include_string(self):
        if self.include_paths:
            return "\n".join([f'include "{p}"' for p in self.include_paths])
        return None

    def get_lib_string(self):
        if self.lib_path_and_names:
            return "\n".join([f'include "{path}" section={libname}' for path, libname in self.lib_path_and_names])
        return None

    def generate_netlist(self, mosfet_model, length, vbs):
        include_str = self.get_include_string()
        lib_str = self.get_lib_string()
        netlist = [
            "simulator lang=spectre",
            "// Lookup Table Generation",
        ]
        if include_str:
            netlist.append(include_str)
        if lib_str:
            netlist.append(lib_str)
        parameters = " ".join(f"{key}={value}" for key, value in self.device_parameters.items())
        netlist.extend([
            "VGS (ng 0) vsource dc=0",
            f"VBS (nb 0) vsource dc={vbs}",
            "VDS (nd 0) vsource dc=0",
            f"{self.mos_spice_symbols[0]} (nd ng 0 nb) {mosfet_model} l={length} {parameters}",
        ])
        if self.raw_spice:
            netlist.extend(self.raw_spice)
        return netlist
