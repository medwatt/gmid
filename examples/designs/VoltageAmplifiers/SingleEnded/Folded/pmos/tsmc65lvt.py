from design import Circuit, run
from mosplot.optimizer import Corner, Knob, Spec

# EDIT: PDK lookup tables and device names.
LUT_DIR = "/home/medwatt/coding/gmid_lookup"
NMOS, PMOS = "nch_lvt", "pch_lvt"
LMIN, LMAX = 130e-9, 2e-6

# EDIT: nominal operating conditions for this PDK.
COND = dict(vdd=1.2, vin_cm=0.6, vout_dc=0.6, cout=5e-12)

# EDIT: knob bounds. Names must match the circuit's KNOBS.
PARAMETERS = [
    Knob("M1a_GMID", (10.0, 20.0)),
    Knob("M2a_GMID", (10.0, 20.0)),
    Knob("M3a_GMID", (10.0, 20.0)),
    Knob("M4a_GMID", (10.0, 20.0)),
    Knob("M5a_GMID", (10.0, 20.0)),
    Knob("M6_GMID",  (10.0, 20.0)),
    Knob("M1a_L", (LMIN, LMAX)),
    Knob("M2a_L", (LMIN, LMAX)),
    Knob("M3a_L", (LMIN, LMAX)),
    Knob("M4a_L", (LMIN, LMAX)),
    Knob("M5a_L", (LMIN, LMAX)),
    Knob("M6_L",  (LMIN, LMAX)),
    Knob("M1a_ID",             (10e-6, 50e-6)),
    Knob("M5a_ID_over_M1a_ID", (1.0, 4.0)),
    Knob("M2a_VDSAT_MARGIN",   (0.05, 0.15), recorner_bound=(-0.5, 0.5)),
    Knob("M5a_VDSAT_MARGIN",   (0.05, 0.15), recorner_bound=(-0.5, 0.5)),
]

# EDIT: targets. Keys must be returned by Circuit.specs().
TARGET_SPECS = {
    "GBW":          Spec(5e6,     "max", 1.0),
    "AC Gain (dB)": Spec(50.0,    "max", 1.0),
    "DC CMR (dB)":  Spec(50.0,    "max", 0.5),
    "PM":           Spec(70.0,    "max", 1.0),
    "Area":         Spec(100e-12, "min", 0.3),
    "Itotal":       Spec(50e-6,  "min", 1.0),
    "Output_Swing": Spec(0.4 * COND["vdd"], "max", 1.0),
    "VOUT_DC":      Spec(COND["vout_dc"], "eq", 1.0),
    "Margin min":   Spec(0.02, "max", 1.0, scale=0.05),
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
