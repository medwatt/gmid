class MosfetNetlistGenerator:
    def __init__(self, model_names, width, mos_spice_symbols, include_paths, lib_path_and_names, raw_spice):
        self.model_names = model_names
        self.width = width
        self.raw_spice = raw_spice
        self.mos_spice_symbols = mos_spice_symbols
        self.include_paths = include_paths
        self.lib_path_and_names = lib_path_and_names

    def get_include_string(self):
        if self.include_paths:
            return "\n".join([f".include '{p}'" for p in self.include_paths])
        return None

    def get_lib_string(self):
        if self.lib_path_and_names:
            return "\n".join([f".lib '{path}' {libname}" for path, libname in self.lib_path_and_names])
        return None

    def generate_netlist(self, mosfet_model, length, vsb):
        include_str = self.get_include_string()
        lib_str = self.get_lib_string()
        transistor_type = self.model_names[mosfet_model]
        r = -1 if transistor_type == "nmos" else 1

        netlist = ["* Lookup Table Generation *"]

        if include_str:
            netlist.append(include_str)

        if lib_str:
            netlist.append(lib_str)

        netlist.extend([
            "VGS NG 0 DC=0",
            f"VBS NB 0 DC={vsb*r}",
            "VDS ND 0 DC=0",
            f"{self.mos_spice_symbols[0]} ND NG 0 NB {mosfet_model} l={length} w={self.width}",
        ])

        if self.raw_spice:
            netlist.extend(self.raw_spice)

        return netlist
