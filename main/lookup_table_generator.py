# -----------------------------------------------------------------------------#
# Author: Mohamed Watfa
# URL: https://github.com/medwatt/
# -----------------------------------------------------------------------------#
import os
import sys
import pickle
import tempfile
import subprocess

import numpy as np

sys.path.append('../parsers')

from parsers.ngspice_parser import NgspiceRawFileReader
from parsers.hspice_parser import import_export

################################################################################

NGSPICE_PATH = "ngspice"
HSPICE_PATH = "hspice"

def range_to_arr(r):
    start, stop, step = r
    return np.arange(start, stop + step, step)

class LookupTableGenerator:
    def __init__(
        self,
        vgs=(0, 1, 0.01),
        vds=(0, 1, 0.01),
        vsb=(0, 1, 0.1),
        width=10e-6,
        lengths=[],
        simulator="ngspice",
        temp=27,
        model_paths=[],
        model_names={"nmos": "NMOS_VTH", "pmos": "PMOS_VTH"},
        description="gmid lookup table",
        raw_spice="",
    ):
        self.vgs = np.array(vgs)
        self.vds = np.array(vds)
        self.vsb = np.array(vsb)
        self.width = width
        self.lengths = np.array(lengths)
        self.simulator = simulator
        self.temp = temp
        self.model_paths = model_paths
        self.model_names = model_names
        self.description = description
        self.raw_spice = raw_spice
        self.lookup_table = {}
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
            "vdsat"
        ]

    ################################################################################

    def __make_tmp_files(self):
        self.input_file_path = tempfile.NamedTemporaryFile(delete=False).name
        self.log_file_path = tempfile.NamedTemporaryFile(delete=False).name
        self.output_file_path = tempfile.NamedTemporaryFile(delete=False).name

    def __remove_tmp_files(self):
        os.remove(self.input_file_path)
        os.remove(self.log_file_path)
        os.remove(self.output_file_path)

    ################################################################################
    #                                   NGSPICE                                    #
    ################################################################################
    def __ngspice_simulator_setup(self, mos):
        save_internal_parameters = "\n".join([f"save @m1[{p}]" for p in self.parameter_names])

        vgs_start, vgs_stop, vgs_step = self.vgs
        vds_start, vds_stop, vds_step = self.vds
        analysis_string = f"dc VDS {vds_start} {vds_stop} {vds_step} VGS {vgs_start} {vgs_stop} {vgs_step}"

        simulator = [
            f".options TEMP = {self.temp}",
            f".options TNOM = {self.temp}",
            ".control",
            save_internal_parameters,
            analysis_string,
            f"write {self.output_file_path} all",
            ".endc",
            ".end",
        ]
        return simulator

    def __run_ngspice(self, circuit):
        with open(self.input_file_path, "w") as file:
            file.write("\n".join(circuit))

        ngspice_command = f"{NGSPICE_PATH} -b -o {self.log_file_path} {self.input_file_path}"
        subprocess.run(ngspice_command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    def __parse_ngspice_output(self):
        ars, _ = NgspiceRawFileReader().read_file(self.output_file_path)
        return ars

    def __save_ngspice_parameters(self, analysis, mos, length, vsb):
        parameter_names = {
            "id": "i(@m1[id])",
            "vth": "v(@m1[vth])",
            "gm": "@m1[gm]",
            "gmbs": "@m1[gmbs]",
            "gds": "@m1[gds]",
            "cgg": "@m1[cgg]",
            "cgs": "@m1[cgs]",
            "cbg": "@m1[cbg]",
            "cgd": "@m1[cgd]",
            "cdd": "@m1[cdd]",
            "vdsat": "v(@m1[vdsat])",
        }

        column_names = analysis[0].dtype.names
        data = analysis[0]

        for p in self.parameter_names:
            col_name = parameter_names[p]
            if col_name in column_names:
                res = np.array(data[col_name])
                self.lookup_table[self.identifier][p][length][vsb] = res.reshape(self.n_vgs, self.n_vds)

    ################################################################################
    #                                    HSPICE                                    #
    ################################################################################
    def __hspice_simulator_setup(self, mos):
        save_internal_parameters = [
            ".probe DC m_id   = par('-i(m1)')",
            ".probe DC m_vt   = par('vth(m1)')",
            ".probe DC m_gm   = par('gmo(m1)')",
            ".probe DC m_gmb  = par('gmbso(m1)')",
            ".probe DC m_gds  = par('gdso(m1)')",
            ".probe DC m_cgg  = par('cggbo(m1)')",
            ".probe DC m_cgs  = par('-cgsbo(m1)')",
            ".probe DC m_cbs  = par('-cbgbo(m1)')",
            ".probe DC m_cgd  = par('-cgdbo(m1)')",
            ".probe DC m_cgb  = par('cggbo(m1)-(-cgsbo(m1))-(-cgdbo(m1))')",
            ".probe DC m_cdd  = par('cddbo(m1)')",
            ".probe DC m_css  = par('-cgsbo(m1)-cbsbo(m1)')",
            ".probe DC m_vdsat = par('vdsat(m1)')",
        ]

        vgs_start, vgs_stop, vgs_step = self.vgs
        vds_start, vds_stop, vds_step = self.vds
        analysis_string = f".dc VGS {vgs_start} {vgs_stop} {vgs_step} VDS {vds_start} {vds_stop} {vds_step}"

        simulator = [
            f".TEMP = {self.temp}",
            ".options dccap brief accurate",
            ".option POST=2",
            "\n".join(save_internal_parameters),
            analysis_string,
            ".end",
        ]
        return simulator

    def __run_hspice(self, circuit):
        with open(self.input_file_path, "w") as file:
            file.write("\n".join(circuit))

        hspice_command = f"{HSPICE_PATH} -i {self.input_file_path} -o {tempfile.gettempdir()}"
        subprocess.run(hspice_command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    def __parse_hspice_output(self):
        import_export(self.input_file_path + ".sw0", "pickle")
        with open(self.input_file_path + "_sw0.pickle", 'rb') as file:
            loaded_data = pickle.load(file)
        return loaded_data

    def __save_hspice_parameters(self, analysis, mos, length, vsb):
        hspice_parameter_names = {
                "id": "i_m1",
                "vth": "vth_m1",
                "gm": "gmo_m1",
                "gmbs": "gmbso_m1",
                "gds": "gdso_m1",
                "cgg": "cggbo_m1",
                "cgs": "cgsbo_m1",
                "cbg": "cbgbo_m1",
                "cgd": "cgdbo_m1",
                "cdd": "cddbo_m1",
                "vdsat": "vdsat_m1",
        }

        for p in self.parameter_names:
            col_name = hspice_parameter_names[p]
            if col_name in analysis.keys():
                res = np.array(analysis[col_name]).T
                self.lookup_table[self.identifier][p][length][vsb] = res

    ################################################################################
    #                                Shared Methods                                #
    ################################################################################
    def __initalize(self):
        self.n_lengths = len(self.lengths)
        self.n_vsb = round((self.vsb[1] - self.vsb[0]) / self.vsb[2]) + 1
        self.n_vds = round((self.vds[1] - self.vds[0]) / self.vds[2]) + 1
        self.n_vgs = round((self.vgs[1] - self.vgs[0]) / self.vgs[2]) + 1

        self.lookup_table[self.identifier] = {}
        for p in self.parameter_names:
            self.lookup_table[self.identifier][p] = np.zeros(shape=(self.n_lengths, self.n_vsb, self.n_vgs, self.n_vds))

        # choose right simulator
        if self.simulator == "ngspice":
            self.simulator_setup = self.__ngspice_simulator_setup
            self.run = self.__run_ngspice
            self.parse = self.__parse_ngspice_output
            self.save = self.__save_ngspice_parameters
        elif self.simulator == "hspice":
            self.simulator_setup = self.__hspice_simulator_setup
            self.run = self.__run_hspice
            self.parse = self.__parse_hspice_output
            self.save = self.__save_hspice_parameters

    def __generate_netlist(self, length, vsb):
        if self.model_paths:
            include_string = "\n".join([f".include '{path}'" for path in self.model_paths])
        else:
            include_string = ""

        circuit = [
            "* Lookup Table Generation *",
            include_string,
            "VGS NG 0 DC=0",
            f"VBS NB 0 DC={-vsb * self.r}",
            "VDS ND 0 DC=0",
            f"M1 ND NG 0 NB {self.model_names[self.identifier]} l={length} w={self.width}",
            self.raw_spice,
        ]
        return circuit

    def __generate_loopkup_table(self, mos):
        self.__initalize()
        for idx, length in enumerate(self.lengths):
            print(f"-- length={length}")
            for idy, vsb in enumerate(np.linspace(self.vsb[0], self.vsb[1], self.n_vsb)):
                circuit = self.__generate_netlist(length, vsb)
                simulator = self.simulator_setup(mos)
                circuit.extend(simulator)
                self.run(circuit)
                analysis = self.parse()
                self.save(analysis, mos, idx, idy)

    def __save_to_dictionary(self):
        self.lookup_table["description"] = self.description
        self.lookup_table["parameter_names"] = self.parameter_names
        self.lookup_table["w"] = self.width
        self.lookup_table["l"] = self.lengths

        if "nmos" in self.model_names:
            self.lookup_table["nmos"]["vgs"] = range_to_arr(self.vgs)
            self.lookup_table["nmos"]["vds"] = range_to_arr(self.vds)
            self.lookup_table["nmos"]["vsb"] = range_to_arr(self.vsb)
            self.lookup_table["nmos"]["model_name"] = self.model_names["nmos"]

        if "pmos" in self.model_names:
            self.lookup_table["pmos"]["vgs"] = -range_to_arr(self.vgs)
            self.lookup_table["pmos"]["vds"] = -range_to_arr(self.vds)
            self.lookup_table["pmos"]["vsb"] = -range_to_arr(self.vsb)
            self.lookup_table["pmos"]["model_name"] = self.model_names["pmos"]

    def __print_netlist(self):
        self.r = 1
        self.identifier = "nmos"
        circuit = self.__generate_netlist(self.lengths[0], 0)
        if self.simulator == "ngspice":
            simulator = self.__ngspice_simulator_setup(self.identifier)
        elif self.simulator == "hspice":
            simulator = self.__hspice_simulator_setup(self.identifier)
        circuit.extend(simulator)
        print("---------------------------------------------------")
        print("----- This is the netlist that gets simulated -----")
        print("---------------------------------------------------")
        print("\n".join(circuit))
        print("---------------------------------------------------")
        print("")

    ################################################################################

    def build(self, filepath):
        self.__make_tmp_files()
        self.__print_netlist()

        if "nmos" in self.model_names:
            print("Generating lookup table for NMOS")
            self.r = 1
            self.identifier = "nmos"
            self.__generate_loopkup_table(self.identifier)

        if "pmos" in self.model_names:
            print("Generating lookup table for PMOS")
            self.r = -1
            self.identifier = "pmos"
            self.__generate_loopkup_table(self.identifier)

        # Save results to file
        print("Saving to file")
        self.__save_to_dictionary()
        np.save(filepath, self.lookup_table, allow_pickle=True)

        # Remove tmp files
        self.__remove_tmp_files()
        print("Done")
