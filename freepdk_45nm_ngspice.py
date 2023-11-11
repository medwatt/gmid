from mosplot import LookupTableGenerator

obj = LookupTableGenerator(
    description="freepdk 45nm ngspice",
    simulator="ngspice",
    model_paths=[
        "/home/medwatt/coding/python/gmid/models/NMOS_VTH.lib",
        "/home/medwatt/coding/python/gmid/models/PMOS_VTH.lib",
        ],
    model_names={
        "nmos": "NMOS_VTH",
        "pmos": "PMOS_VTH",
    },
    vsb=(0, 1.0, 0.1),
    vgs=(0, 1.0, 0.01),
    vds=(0, 1.0, 0.01),
    width=10e-6,
    lengths=[45e-9, 100e-9, 200e-9, 400e-9, 800e-9, 1.6e-6, 3.2e-6, 6.4e-6],
)
obj.build("/home/medwatt/coding/python/gmid/lookup_tables/freepdk_45nm_ngspice.npy")
