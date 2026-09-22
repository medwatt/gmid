from design import Circuit, run
from mosplot.optimizer import Corner, Knob, Spec

# EDIT: PDK lookup tables and device names.
LUT_DIR = "/home/medwatt/coding/gmid_lookup"
NMOS, PMOS = "sg13_hv_nmos", "sg13_hv_pmos"

# EDIT: nominal operating conditions for this PDK.
COND = dict(vdd=3.3)

# EDIT: knob bounds. Names must match the circuit's KNOBS.
PARAMETERS = [
    Knob("M2_GMID", (8.0, 18.0)),
    Knob("M3_GMID", (8.0, 18.0)),
    Knob("L_n", (500e-9, 5e-6)),
    Knob("L_p", (500e-9, 5e-6)),
    Knob("K", (4.0, 8.0)),
    Knob("IREF", (10e-6, 60e-6)),
]

# EDIT: targets. Keys must be returned by Circuit.specs().
TARGET_SPECS = {
    "gm_ref": Spec(100e-6, "max", 3.0),
    "gmR": Spec(1.0, "eq", 20.0),
    "VDS_margin": Spec(0.15, "max", 10.0),
    "Area": Spec(200e-12, "min", 0.3),
    "Itotal": Spec(80e-6, "min", 0.5),
}

# EDIT: one Corner per process LUT (voltage corners reuse a LUT with overridden COND).
CORNERS = [
    Corner("tt", f"{LUT_DIR}/ihp_130_hv_spectre.npz", NMOS, PMOS, conditions=COND),
]

if __name__ == "__main__":
    run(
        workers=4,
        circuit=Circuit,
        parameters=PARAMETERS,
        target_specs=TARGET_SPECS,
        corners=CORNERS,
        output_module="./design.scs",
        maxiter=60,
        seed=1,
        n_restarts=1,
    )
