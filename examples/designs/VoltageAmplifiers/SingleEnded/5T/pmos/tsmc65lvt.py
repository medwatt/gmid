from design import Circuit, run
from mosplot.optimizer import Corner, Knob, Spec

# EDIT: PDK lookup tables and device names.
LUT_DIR = "/home/medwatt/coding/gmid_lookup"
NMOS, PMOS = "nch_lvt", "pch_lvt"
LMIN, LMAX = 100e-9, 5e-6

# EDIT: nominal operating conditions for this PDK.
COND = dict(vdd=1.2, vin_cm=0.4, vout_dc=0.6, cout=5e-12)

# EDIT: knob bounds. Names must match the circuit's KNOBS.
PARAMETERS = [
    Knob("M1a_GMID", (10, 20)),
    Knob("M2a_GMID", (10, 20)),
    Knob("M3_GMID", (8, 20)),
    Knob("M1a_L", (LMIN, LMAX)),
    Knob("M2a_L", (LMIN, LMAX)),
    Knob("M3_L", (LMIN, LMAX)),
    Knob("M1a_ID", (5e-06, 2e-05)),
]

# EDIT: targets. Keys must be returned by Circuit.specs().
TARGET_SPECS = {
    "GBW": Spec(5e6, "max", 5.0),
    "AC Gain (dB)": Spec(30.0, "max", 0.5),
    "DC CMR (dB)": Spec(40.0, "max", 0.5),
    "PM": Spec(90.0, "max", 0.1),
    "Area": Spec(10e-12, "min", 1.0),
    "Itotal": Spec(20e-6, "min", 2.0),
    "Output_Swing": Spec(0.5 * COND["vdd"], "max", 1.0),
    "VOUT_DC": Spec(COND["vout_dc"], "eq", 2.0),
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
        maxiter=200,
        seed=1,
        n_restarts=1,
    )
