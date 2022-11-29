#-----------------------------------------------------------------------------#
# Author: Mohamed Watfa
# URL: https://github.com/medwatt/
#-----------------------------------------------------------------------------#
import numpy as np
from PySpice.Spice.Netlist import Circuit
from PySpice.Spice.Library import SpiceLibrary


def range_to_arr(r):
    start, stop, step = r
    return np.arange(start, stop + step, step)


class LookupTableGenerator:
    def __init__(
        self,
        vgs=(0, 1, 0.01),
        vds=(0, 1, 0.01),
        vsb=(0, 1, 0.01),
        width=10e-6,
        lengths=(1e-6),
        temp=27,
        models_path="./models",
        model_names={"nmos": "NMOS_VTH", "pmos": "PMOS_VTH"},
        description="gmid lookup table",
    ):
        self.vgs = np.array(vgs)
        self.vds = np.array(vds)
        self.vsb = np.array(vsb)
        self.width = width
        self.lengths = np.array(lengths)
        self.temp = temp
        self.spice_library = SpiceLibrary(models_path)
        self.model_names = model_names
        self.description = description
        self.parameter_names = [
            "id",
            "vth",
            "gm",
            "gmbs",
            "gds",
            "cgg",
            "cgs",
            "cbg",
            "cgd",
            "cdd",
        ]
        self.lookup_table = {}

    def __initalize(self):
        self.n_lengths = len(self.lengths)
        self.n_vsb = round((self.vsb[1] - self.vsb[0]) / self.vsb[2]) + 1
        self.n_vds = round((self.vds[1] - self.vds[0]) / self.vds[2]) + 1
        self.n_vgs = round((self.vgs[1] - self.vgs[0]) / self.vgs[2]) + 1
        self.lookup_table[self.identifier] = {}
        for p in self.parameter_names:
            self.lookup_table[self.identifier][p] = np.zeros(
                shape=(self.n_lengths, self.n_vsb, self.n_vds, self.n_vgs)
            )

    def __build_netlist(self, length, vsb):
        circuit = Circuit("gm/id simulation")
        circuit.V("GS", "NG", circuit.gnd, 0)
        circuit.V("BS", "NB", circuit.gnd, -vsb * self.r)
        circuit.V("DS", "ND", circuit.gnd, 0)
        circuit.include(self.spice_library[self.model_names[self.identifier]])
        circuit.M(
            self.identifier,
            "ND",
            "NG",
            circuit.gnd,
            "NB",
            model=self.model_names[self.identifier],
            l=length,
            w=self.width,
        )
        return circuit

    def __simulate_circuit(self, circuit, mos):
        simulator = circuit.simulator(
            temperature=self.temp, nominal_temperature=self.temp
        )
        simulator.save_internal_parameters(
            *[f"@{mos}[{p}]" for p in self.parameter_names]
        )
        analysis = simulator.dc(
            VGS=slice(*(self.vgs * self.r)), VDS=slice(*(self.vds * self.r))
        )
        return analysis

    def __save_parameters(self, analysis, mos, length, vsb):
        for p in self.parameter_names:
            res = analysis[f"@{mos}[{p}]"]
            self.lookup_table[self.identifier][p][length][vsb] = res.reshape(
                self.n_vds, self.n_vgs
            )

    def __generate_loopkup_table(self, mos):
        self.__initalize()
        for idx, length in enumerate(self.lengths):
            for idy, vsb in enumerate(
                np.linspace(self.vsb[0], self.vsb[1], self.n_vsb)
            ):
                circuit = self.__build_netlist(length, vsb)
                analysis = self.__simulate_circuit(circuit, mos)
                self.__save_parameters(analysis, mos, idx, idy)

    def __save_to_dictionary(self):
        self.lookup_table["nmos"]["vgs"] = range_to_arr(self.vgs)
        self.lookup_table["nmos"]["vds"] = range_to_arr(self.vds)
        self.lookup_table["nmos"]["vsb"] = range_to_arr(self.vsb)
        self.lookup_table["nmos"]["w"] = self.width
        self.lookup_table["nmos"]["l"] = self.lengths
        self.lookup_table["nmos"]["parameter_names"] = self.parameter_names
        self.lookup_table["nmos"]["model_name"] = self.model_names["nmos"]
        self.lookup_table["pmos"]["vgs"] = -range_to_arr(self.vgs)
        self.lookup_table["pmos"]["vds"] = -range_to_arr(self.vds)
        self.lookup_table["pmos"]["vsb"] = -range_to_arr(self.vsb)
        self.lookup_table["pmos"]["w"] = self.width
        self.lookup_table["pmos"]["l"] = self.lengths
        self.lookup_table["pmos"]["parameter_names"] = self.parameter_names
        self.lookup_table["pmos"]["model_name"] = self.model_names["pmos"]
        self.lookup_table["description"] = self.description

    def build(self, filepath):
        # NMOS
        print("Generating lookup table for NMOS")
        self.r = 1
        self.identifier = "nmos"
        self.__generate_loopkup_table("m" + self.identifier)
        # PMOS
        print("Generating lookup table for PMOS")
        self.r = -1
        self.identifier = "pmos"
        self.__generate_loopkup_table("m" + self.identifier)
        # save results to dictionary and to file
        print("Saving to file")
        self.__save_to_dictionary()
        np.save(filepath, self.lookup_table, allow_pickle=True)


if __name__ == "__main__":

    # example
    obj = LookupTableGenerator (
        vgs=(0, 1, 0.01),
        vds=(0, 1, 0.05),
        vsb=(0, 1, 0.1),
        width=10e-6,
        lengths=[50e-9, 100e-9, 200e-9, 400e-9, 800e-9, 1.6e-6, 3.2e-6, 6.4e-6],
        models_path="./models",
        model_names={
            "nmos": "NMOS_VTH",
            "pmos": "PMOS_VTH"},
        description="freepdk 45nm"
        )
    obj.build("./freepdk45_loopup_table.npy")
