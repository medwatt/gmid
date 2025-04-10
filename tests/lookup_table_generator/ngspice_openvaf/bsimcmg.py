from mosplot.lookup_table_generator import LookupTableGenerator, TransistorSweep
from mosplot.lookup_table_generator.simulators import NgspiceSimulator

ngspice = NgspiceSimulator(
    include_paths=[
        "./Modelcards/modelcard.nmos",
        "./Modelcards/modelcard.pmos",
    ],
    osdi_paths = [
        "./bsimcmg.osdi"
    ],
    mos_spice_symbols = ("n1", "n1"),
    device_parameters = {
        "NF": 10,
    }
)

nmos_sweep = TransistorSweep(
    vgs=(0, 1.0, 0.01),
    vds=(0, 1.0, 0.01),
    vbs=(-1.0, 1.0, 0.1),
    length=[100e-9, 200e-9],
)

pmos_sweep = TransistorSweep(
    vgs=(0, -1.0, -0.01),
    vds=(0, -1.0, -0.01),
    vbs=(-1.0, 1.0, 0.1),
    length=[100e-9, 200e-9],
)

obj = LookupTableGenerator(
    description="BSIMCMG test",
    simulator=ngspice,
    model_sweeps={
        "BSIMCMG_osdi_N": nmos_sweep,
        "BSIMCMG_osdi_P": pmos_sweep,
    },
)

# obj.op_simulation()
obj.build("./bsimcmg")
