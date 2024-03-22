# -----------------------------------------------------------------------------#
# Author: Mohamed Watfa
# URL: https://github.com/medwatt/
# -----------------------------------------------------------------------------#

import os
import shutil
import pickle
import tempfile
import subprocess

import numpy as np

from ..parsers.ngspice_parser import NgspiceRawFileReader
from ..parsers.hspice_parser import import_export

################################################################################

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
        lengths=[500e-9, 600e-9],
        simulator="ngspice",
        simulator_path=None,
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
        self.simulator_path = simulator_path
        self.temp = temp
        self.model_paths = model_paths
        self.model_names = model_names
        self.description = description
        self.raw_spice = raw_spice
        self.lookup_table = {}

        self.validate_paths()
        self.setup_simulator()

    ################################################################################
    def validate_paths(self):
        path_list = self.model_paths[::]
        if self.simulator_path is not None:
            path_list.append(self.simulator_path)
        for path in path_list:
            if not os.path.exists(path):
                raise FileNotFoundError(f"The path {path} does not exist.")

    def setup_simulator(self):
        def check_for_binary(binary_name):
            if shutil.which(binary_name) is None:
                raise ValueError(f"The binary '{binary_name}' is not accessible.")
        if self.simulator_path is None:
            check_for_binary(self.simulator)
            self.simulator_path = self.simulator

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
    def __ngspice_parameters(self):
        self.parameter_table = {
            # parameter name : [name recognized by simulator, name used in the output file],
            "id"   : ["save i(vds)"    , "i(i_vds)"],
            "vth"  : ["save @m1[vth]"  , "v(@m1[vth])"],
            "vdsat": ["save @m1[vdsat]", "v(@m1[vdsat])"],
            "gm"   : ["save @m1[gm]"   , "@m1[gm]"],
            "gmbs" : ["save @m1[gmbs]" , "@m1[gmbs]"],
            "gds"  : ["save @m1[gds]"  , "@m1[gds]"],
            "cgg"  : ["save @m1[cgg]"  , "@m1[cgg]"],
            "cgs"  : ["save @m1[cgs]"  , "@m1[cgs]"],
            "cbg"  : ["save @m1[cbg]"  , "@m1[cbg]"],
            "cgd"  : ["save @m1[cgd]"  , "@m1[cgd]"],
            "cdd"  : ["save @m1[cdd]"  , "@m1[cdd]"],
        }
        self.save_internal_parameters = "\n".join([values[0] for values in self.parameter_table.values()])

    def __ngspice_simulator_setup(self):
        vgs_start, vgs_stop, vgs_step = self.vgs * self.r
        vds_start, vds_stop, vds_step = self.vds * self.r
        analysis_string = f"dc VDS {vds_start} {vds_stop} {vds_step} VGS {vgs_start} {vgs_stop} {vgs_step}"

        simulator = [
            f".options TEMP = {self.temp}",
            f".options TNOM = {self.temp}",
            ".control",
            self.save_internal_parameters,
            analysis_string,
            "let i_vds = abs(i(vds))",
            f"write {self.output_file_path} all",
            ".endc",
            ".end",
        ]
        return simulator

    def __run_ngspice(self, circuit):
        with open(self.input_file_path, "w") as file:
            file.write("\n".join(circuit))
        ngspice_command = f"{self.simulator_path} -b -o {self.log_file_path} {self.input_file_path}"
        subprocess.run(ngspice_command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    def __parse_ngspice_output(self):
        ars, _ = NgspiceRawFileReader().read_file(self.output_file_path)
        return ars

    def __save_ngspice_parameters(self, analysis, mos, length, vsb):
        column_names = analysis[0].dtype.names
        data = analysis[0]

        for p in self.parameter_table.keys():
            col_name = self.parameter_table[p][1]
            if col_name in column_names:
                res = np.array(data[col_name])
                self.lookup_table[self.identifier][p][length][vsb] = res.reshape(self.n_vgs, self.n_vds)

    ################################################################################
    #                                    HSPICE                                    #
    ################################################################################
    def __hspice_parameters(self):
        self.parameter_table = {
            # parameter name : [name recognized by simulator, name used in the output file],
            "id"   : [".probe DC m_id    = par('abs(i(vds))')"                        , "m_id"],
            "vth"  : [".probe DC m_vth   = par('vth(m1)')"                            , "m_vth"],
            "vdsat": [".probe DC m_vdsat = par('vdsat(m1)')"                          , "m_vdsat"],
            "gm"   : [".probe DC m_gm    = par('gmo(m1)')"                            , "m_gm"],
            "gmbs" : [".probe DC m_gmb   = par('gmbso(m1)')"                          , "m_gmb"],
            "gds"  : [".probe DC m_gds   = par('gdso(m1)')"                           , "m_gds"],
            "cgg"  : [".probe DC m_cgg   = par('cggbo(m1)')"                          , "m_cgg"],
            "cgs"  : [".probe DC m_cgs   = par('-cgsbo(m1)')"                         , "m_cgs"],
            "cgd"  : [".probe DC m_cgd   = par('-cgdbo(m1)')"                         , "m_cgd"],
            "cgb"  : [".probe DC m_cgb   = par('cggbo(m1)-(-cgsbo(m1))-(-cgdbo(m1))')", "m_cgb"],
            "cdd"  : [".probe DC m_cdd   = par('cddbo(m1)')"                          , "m_cdd"],
            "css"  : [".probe DC m_css   = par('-cgsbo(m1)-cbsbo(m1)')"               , "m_css"],
        }
        self.save_internal_parameters = "\n".join([values[0] for values in self.parameter_table.values()])

    def __hspice_simulator_setup(self):
        vgs_start, vgs_stop, vgs_step = self.vgs * self.r
        vds_start, vds_stop, vds_step = self.vds * self.r
        analysis_string = f".dc VGS {vgs_start} {vgs_stop} {vgs_step} VDS {vds_start} {vds_stop} {vds_step}"

        simulator = [
            f".TEMP = {self.temp}",
            ".options probe dccap brief accurate",
            ".option POST=2",
            self.save_internal_parameters,
            analysis_string,
            ".end",
        ]
        return simulator

    def __run_hspice(self, circuit):
        with open(self.input_file_path, "w") as file:
            file.write("\n".join(circuit))
        hspice_command = f"{self.simulator_path} -i {self.input_file_path} -o {tempfile.gettempdir()}"
        subprocess.run(hspice_command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    def __parse_hspice_output(self):
        import_export(self.input_file_path + ".sw0", "pickle")
        with open(self.input_file_path + "_sw0.pickle", 'rb') as file:
            loaded_data = pickle.load(file)
        return loaded_data

    def __save_hspice_parameters(self, analysis, mos, length, vsb):
        for p in self.parameter_table.keys():
            col_name = self.parameter_table[p][1]
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
        for p in self.parameter_table:
            self.lookup_table[self.identifier][p] = np.zeros(shape=(self.n_lengths, self.n_vsb, self.n_vgs, self.n_vds))

        # choose right simulator
        if self.simulator == "ngspice":
            self.__ngspice_parameters()
            self.simulator_setup = self.__ngspice_simulator_setup
            self.run = self.__run_ngspice
            self.parse = self.__parse_ngspice_output
            self.save = self.__save_ngspice_parameters
        elif self.simulator == "hspice":
            self.__hspice_parameters()
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
                simulator = self.simulator_setup()
                circuit.extend(simulator)
                self.run(circuit)
                analysis = self.parse()
                self.save(analysis, mos, idx, idy)

    def __save_to_dictionary(self):
        self.lookup_table["description"] = self.description
        self.lookup_table["parameter_names"] = list(self.parameter_table.keys())
        self.lookup_table["width"] = self.width
        self.lookup_table["lengths"] = self.lengths

        if "nmos" in self.model_names:
            self.lookup_table["nmos"]["vgs"] = range_to_arr(self.vgs)
            self.lookup_table["nmos"]["vds"] = range_to_arr(self.vds)
            self.lookup_table["nmos"]["vsb"] = range_to_arr(self.vsb)
            self.lookup_table["nmos"]["model_name"] = self.model_names["nmos"]
            self.lookup_table["nmos"]["width"] = self.width
            self.lookup_table["nmos"]["lengths"] = self.lengths
            self.lookup_table["nmos"]["parameter_names"] = list(self.parameter_table.keys())

        if "pmos" in self.model_names:
            self.lookup_table["pmos"]["vgs"] = -range_to_arr(self.vgs)
            self.lookup_table["pmos"]["vds"] = -range_to_arr(self.vds)
            self.lookup_table["pmos"]["vsb"] = -range_to_arr(self.vsb)
            self.lookup_table["pmos"]["model_name"] = self.model_names["pmos"]
            self.lookup_table["pmos"]["width"] = self.width
            self.lookup_table["pmos"]["lengths"] = self.lengths
            self.lookup_table["pmos"]["parameter_names"] = list(self.parameter_table.keys())

    def __print_netlist(self):
        self.r = 1
        self.identifier = "nmos"
        circuit = self.__generate_netlist(self.lengths[0], 0)
        if self.simulator == "ngspice":
            self.__ngspice_parameters()
            simulator = self.__ngspice_simulator_setup()
        elif self.simulator == "hspice":
            self.__hspice_parameters()
            simulator = self.__hspice_simulator_setup()
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
