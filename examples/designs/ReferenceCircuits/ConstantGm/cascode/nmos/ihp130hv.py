from design import Circuit, run
from mosplot.optimizer import Corner, Knob, Spec

LUT_DIR = "/home/medwatt/coding/gmid_lookup"
NMOS, PMOS = "sg13_hv_nmos", "sg13_hv_pmos"
COND = dict(vdd=3.3)

PARAMETERS = [
    Knob("M2_GMID", (8.0, 18.0)),
    Knob("M4_GMID", (6.0, 16.0)),
    Knob("M5_GMID", (6.0, 16.0)),
    Knob("M7_GMID", (8.0, 18.0)),
    Knob("L_n", (500e-9, 5e-6)),
    Knob("L_nc", (400e-9, 5e-6)),
    Knob("L_pc", (400e-9, 5e-6)),
    Knob("L_pm", (500e-9, 5e-6)),
    Knob("K", (4.0, 8.0)),
    Knob("IREF", (10e-6, 60e-6)),
]

TARGET_SPECS = {
    "gm_ref": Spec(100e-6, "max", 3.0),  # want gm_ref >= 100 uS
    "gmR": Spec(1.0, "eq", 20.0),  # constant-gm figure of merit
    "VDS_margin": Spec(0.15, "max", 10.0),  # keep every cascode >=150 mV into saturation
    "Area": Spec(500e-12, "min", 0.3),
    "Itotal": Spec(80e-6, "min", 0.5),
}

CORNERS = [
    Corner("tt", f"{LUT_DIR}/ihp_130_hv_spectre.npz", NMOS, PMOS, conditions=COND),
    Corner("ss", f"{LUT_DIR}/ihp_130_hv_spectre_ss.npz", NMOS, PMOS, conditions=COND),
    Corner("ff", f"{LUT_DIR}/ihp_130_hv_spectre_ff.npz", NMOS, PMOS, conditions=COND),
]

if __name__ == "__main__":
    run(
        workers=4,
        circuit=Circuit,
        parameters=PARAMETERS,
        target_specs=TARGET_SPECS,
        corners=CORNERS,
        output_module="./design.scs",
        maxiter=80,
        seed=1,
        n_restarts=1,
    )
