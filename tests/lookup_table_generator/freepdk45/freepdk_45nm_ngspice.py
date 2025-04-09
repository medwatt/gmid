from mosplot.lookup_table_generator import LookupTableGenerator, TransistorSweep
from mosplot.lookup_table_generator.simulators import NgspiceSimulator, HspiceSimulator

# Create an ngspice simulator instance.
ngspice = NgspiceSimulator(
    simulator_path="ngspice",
    temperature=25,
    parameters_to_save=["id", "vth", "vdsat", "gm"],
    mos_spice_symbols=("m1", "m1"),
    include_paths=[
        "./NMOS_VTH.inc",
        "./PMOS_VTH.inc"
    ],
)

# Create an hspice simulator instance.
# hspice = HspiceSimulator(
#     simulator_path="hspice",
#     temperature=25,
#     parameters_to_save=["id", "vth", "vdsat", "gm"],
#     mos_spice_symbols=("m1", "m1"),
#     include_paths=[
#         "./NMOS_VTH.inc",
#         "./PMOS_VTH.inc"
#     ],
# )

# Define separate sweep objects for NMOS transistors.
nmos_sweep = TransistorSweep(
    vgs=(0, 1.0, 0.01),
    vds=(0, 1.0, 0.01),
    vbs=(0, -1.0, -0.1),
    length=[45e-9, 100e-9]
)

# Define separate sweep objects for PMOS transistors.
pmos_sweep = TransistorSweep(
    vgs=(0, -1.0, -0.01),
    vds=(0, -1.0, -0.01),
    vbs=(0, 1.0, 0.1),
    length=[200e-9, 500e-9]
)

# Create a lookup table generator.
obj = LookupTableGenerator(
    description="freepdk 45nm",
    simulator=ngspice,
    model_sweeps={
        "NMOS_VTH": nmos_sweep,
        "PMOS_VTH": pmos_sweep,
    },
    width=10e-6,
)

obj.build("./freepdk45")
# obj.op_simulation()
