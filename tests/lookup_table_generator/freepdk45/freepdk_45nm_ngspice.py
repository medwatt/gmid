from mosplot.lookup_table_generator import LookupTableGenerator, TransistorSweep

# Define separate sweep objects for NMOS.
nmos_sweep = TransistorSweep(
    vgs=(0, 1.0, 0.01),
    vds=(0, 1.0, 0.01),
    vbs=(0, -1.0, -0.1),
    length=[45e-9, 100e-9]
)

# Define separate sweep objects for PMOS.
pmos_sweep = TransistorSweep(
    vgs=(0, -1.0, -0.01),
    vds=(0, -1.0, -0.01),
    vbs=(0, 1.0, 0.1),
    length=[200e-9, 500e-9]
)

obj = LookupTableGenerator(
    description="freepdk 45nm",
    simulator="ngspice",
    include_paths=[
        "./NMOS_VTH.inc",
        "./PMOS_VTH.inc",
        ],
    model_sweeps={
        "NMOS_VTH": nmos_sweep,
        "PMOS_VTH": pmos_sweep,
    },
    width=10e-6,
)

obj.build("./freepdk45")
# obj.op_simulation()
