from mosplot.lookup_table_generator import LookupTableGenerator, TransistorSweep
from mosplot.lookup_table_generator.simulators import NgspiceSimulator

# Create an ngspice simulator instance.
ngspice = NgspiceSimulator(
    include_paths=[
        "./NMOS_VTH.inc",
        "./PMOS_VTH.inc"
    ],
    device_parameters = {
        "w": 10e-6,
    }
)

# Define a sweep object for NMOS transistors.
nmos_sweep = TransistorSweep(
    mos_type = "nmos",
    vgs=(0, 1.0, 0.01),
    vds=(0, 1.0, 0.01),
    vbs=(0, -1.0, -0.1),
    length=[45e-9, 100e-9],
)

# Define a sweep object for PMOS transistors.
pmos_sweep = TransistorSweep(
    mos_type = "pmos",
    vgs=(0, -1.0, -0.01),
    vds=(0, -1.0, -0.01),
    vbs=(0, 1.0, 0.1),
    length=[200e-9, 500e-9],
)

# Create a lookup table generator.
obj = LookupTableGenerator(
    description="freepdk 45nm",
    simulator=ngspice,
    model_sweeps={
        "NMOS_VTH": nmos_sweep,
        "PMOS_VTH": pmos_sweep,
    },
    n_process=1,
)

# obj.op_simulation()
obj.build("./freepdk45")
