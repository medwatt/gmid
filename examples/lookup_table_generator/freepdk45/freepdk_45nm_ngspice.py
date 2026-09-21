import numpy as np

from mosplot.lookup_table_generator import LookupTableGenerator, TransistorSweep
from mosplot.lookup_table_generator.simulators import NgspiceSimulator, HspiceSimulator, SpectreSimulator

SIMULATOR = NgspiceSimulator

DESCRIPTION = "Freepdk 45nm GMID lookup table"

BUILD_DIR = "/home/medwatt/freepdk_45_ngspice"

INCLUDE_PATHS=[
    "NMOS_VTH.inc",
    "PMOS_VTH.inc",
]

NMOS = "NMOS_VTH"
PMOS = "PMOS_VTH"

LENGTHS = np.array([
    0.45, 0.1,  0.15, 0.20, 0.24,
    0.28, 0.33, 0.39, 0.47, 0.56,
    0.68, 0.82, 1.00, 1.20, 1.50,
    1.80, 2.20, 2.70, 3.30, 4.00,
    5.00, 6.20, 7.50, 9.00, 10.00,
]) * 1e-6

VDD = 1.0
VGS_DELTA = 0.01
VDS_DELTA = 0.01
VBS_DELTA = 0.1

def main():

    sim = SIMULATOR(
        include_paths=INCLUDE_PATHS,
        device_parameters={
            "w": 10e-6,
        },
    )

    nmos_sweep = TransistorSweep(
        mos_type="nmos",
        vgs=(0.0, VDD, VGS_DELTA),
        vds=(VDS_DELTA, VDD, VDS_DELTA),
        vbs=(0.0, -VDD, -VBS_DELTA),
        length=LENGTHS,
    )

    pmos_sweep = TransistorSweep(
        mos_type="pmos",
        vgs=(0.0, -VDD, -VGS_DELTA),
        vds=(-VDS_DELTA, -VDD, -VDS_DELTA),
        vbs=(0.0,  VDD,  VBS_DELTA),
        length=LENGTHS,
    )

    obj = LookupTableGenerator(
        description=DESCRIPTION,
        simulator=sim,
        model_sweeps={
            NMOS: nmos_sweep,
            PMOS: pmos_sweep,
        },
        n_process=1,
    )

    obj.build(BUILD_DIR)


if __name__ == "__main__":
    main()
