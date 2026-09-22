from design import Circuit, run
from mosplot.optimizer import Corner, Knob, Spec

# EDIT: PDK lookup tables and device names.
LUT_DIR = "/home/medwatt/coding/gmid_lookup"
NMOS, PMOS = "nch_lvt", "pch_lvt"
LMIN, LMAX = 130e-9, 2e-6

# EDIT: nominal operating conditions for this PDK.
COND = dict(vdd=1.2, vin_cm=0.9, vout_dc=0.6, cout=5e-12)

# EDIT: knob bounds. Names must match the circuit's KNOBS.
PARAMETERS = [
    Knob("M1a_GMID", (10, 20)),
    Knob("M2a_GMID", (10, 20)),
    Knob("M3_GMID", (10, 20)),
    Knob("M4_GMID", (10, 20)),
    Knob("M1a_L", (LMIN, LMAX)),
    Knob("M2a_L", (LMIN, LMAX)),
    Knob("M3_L", (LMIN, LMAX)),
    Knob("M4_L", (LMIN, LMAX)),
    Knob("M5_over_M3", (1, 16)),
    Knob("CC", (1e-12, 2e-12)),
    Knob("Rz", (100, 10000)),
    Knob("M1a_ID", (5e-06, 30e-06)),
]

# EDIT: targets. Keys must be returned by Circuit.specs().
TARGET_SPECS = {
    "GBW": Spec(10e6, "max", 1.0),
    "AC Gain (dB)": Spec(50.0, "max", 1.0),
    "DC CMR (dB)": Spec(50.0, "max", 0.5),
    "PM": Spec(70.0, "max", 2.0),
    "Area": Spec(50e-12, "min", 1.0),
    "Itotal": Spec(50e-6, "min", 1.0),


    "VOUT_Error": Spec(0.02, "min", 20.0),
}

# EDIT: one Corner per process LUT (voltage corners reuse a LUT with overridden COND).
CORNERS = [
    Corner("tt", f"{LUT_DIR}/tsmc_65_lv_spectre_tt.npz", NMOS, PMOS, conditions=COND),
    Corner("ff", f"{LUT_DIR}/tsmc_65_lv_spectre_ff.npz", NMOS, PMOS, conditions=COND),
    Corner("ss", f"{LUT_DIR}/tsmc_65_lv_spectre_ss.npz", NMOS, PMOS, conditions=COND),
]

if __name__ == "__main__":
    run(
        workers=4,
        circuit=Circuit,
        parameters=PARAMETERS,
        target_specs=TARGET_SPECS,
        corners=CORNERS,
        output_module="./design.scs",
        maxiter=250,
        seed=1,
        n_restarts=1,
    )
