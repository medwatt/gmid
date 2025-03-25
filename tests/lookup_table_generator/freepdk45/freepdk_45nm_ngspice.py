from mosplot.lookup_table_generator import LookupTableGenerator

obj = LookupTableGenerator(
    description="freepdk 45nm",
    simulator="ngspice",
    include_paths=[
        "./NMOS_VTH.inc",
        "./PMOS_VTH.inc",
        ],
    model_names={
        "NMOS_VTH": "nmos",
        "PMOS_VTH": "pmos",
    },
    vsb=(0, 1.0, 0.1),
    vgs=(0, 1.0, 0.01),
    vds=(0, 1.0, 0.01),
    width=10e-6,
    lengths=[45e-9, 100e-9],
)

# obj.op_simulation()

obj.build("./freepdk_45nm_ngspice")
