from design import Circuit, run
from mosplot.optimizer import Corner, Knob, Spec

LUT_DIR = "/home/medwatt/coding/gmid_lookup"
NMOS, PMOS = "nch_lvt", "pch_lvt"

# EDIT: nominal operating conditions for this PDK.
COND = dict(vdd=1.2, vin_cm=0.5, vout_dc=0.6, cout=5e-12)
LMIN, LMAX = 130e-9, 2e-6

PARAMETERS = [
    Knob("M1a_GMID", (10, 20)),
    Knob("M2a_GMID", (10, 20)),
    Knob("M3a_GMID", (10, 20)),
    Knob("M4a_GMID", (10, 20)),
    Knob("M5_GMID", (10, 20)),
    Knob("M1a_L", (LMIN, LMAX)),
    Knob("M2a_L", (LMIN, LMAX)),
    Knob("M3a_L", (LMIN, LMAX)),
    Knob("M4a_L", (LMIN, LMAX)),
    Knob("M5_L", (LMIN, LMAX)),
    Knob("M1a_ID", (5e-06, 40e-06)),
    Knob("M4a_VDSAT_MARGIN", (0.05, 0.26), recorner_bound=(-0.5, 0.5)),
    Knob("M1a_VDSAT_MARGIN", (0.05, 0.26), recorner_bound=(-0.5, 0.5)),
]

# EDIT: targets. Keys must be returned by Circuit.specs().
TARGET_SPECS = {
    "GBW": Spec(10e6, "max", 1.0),
    "AC Gain (dB)": Spec(50.0, "max", 1.0),
    "DC CMR (dB)": Spec(50.0, "max", 0.5),
    "PM": Spec(70.0, "max", 2.0),
    "Area": Spec(50e-12, "min", 1.0),
    "Itotal": Spec(20e-6, "min", 1.0),
    "Output_Swing": Spec(0.4 * COND["vdd"], "max", 1.0),
    "VOUT_DC": Spec(COND["vout_dc"], "eq", 2.0),
    "Margin min":   Spec(0.02, "max", 1.0, scale=0.05),   # saturation margin of the devices under the fixed biases, every corner
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
        maxiter=100,
        seed=1,
        n_restarts=1,
    )
