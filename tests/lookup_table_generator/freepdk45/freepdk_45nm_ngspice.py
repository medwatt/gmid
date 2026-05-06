from mosplot.lookup_table_generator import LookupTableGenerator, TransistorSweep
from mosplot.lookup_table_generator.simulators import NgspiceSimulator

def main():
    ngspice = NgspiceSimulator(
        include_paths=[
            "/home/medwatt/git/gmid/tests/lookup_table_generator/freepdk45/NMOS_VTH.inc",
            "/home/medwatt/git/gmid/tests/lookup_table_generator/freepdk45/PMOS_VTH.inc"
        ],
        device_parameters={"w": 10e-6},
    )

    nmos_sweep = TransistorSweep(
        mos_type="nmos",
        vgs=(0, 1.0, 0.01),
        vds=(0, 1.0, 0.01),
        vbs=(0, -1.0, -0.5),
        length=[45e-9, 100e-9],
    )

    pmos_sweep = TransistorSweep(
        mos_type="pmos",
        vgs=(0, -1.0, -0.01),
        vds=(0, -1.0, -0.01),
        vbs=(0, 1.0, 0.5),
        length=[45e-9, 100e-9],
    )

    obj = LookupTableGenerator(
        description="freepdk 45nm",
        simulator=ngspice,
        model_sweeps={"NMOS_VTH": nmos_sweep, "PMOS_VTH": pmos_sweep},
        n_process=1,
    )

    obj.build("/home/medwatt/git/gmid/tests/lookup_table_generator/local_testing/builds/freepdk45_ngspice")

if __name__ == "__main__":
    main()
