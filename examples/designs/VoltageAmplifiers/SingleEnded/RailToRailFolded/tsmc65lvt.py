from design import Circuit, run
from mosplot.optimizer import Corner, Knob, Spec

# EDIT: PDK lookup tables and device names.
LUT_DIR = "/home/medwatt/coding/gmid_lookup"
NMOS, PMOS = "nch_lvt", "pch_lvt"
LMIN, LMAX = 100e-9, 1e-6

# EDIT: nominal operating conditions for this PDK.
COND = dict(vdd=1.2, vin_cm=0.6, vout_dc=0.6, cout=5e-12)

# EDIT: knob bounds. Names must match the circuit's KNOBS.
PARAMETERS = [
    Knob("MP1a_GMID", (5, 25)),
    Knob("MN1a_GMID", (5, 25)),
    Knob("M2a_GMID", (5, 25)),
    Knob("M3a_GMID", (5, 25)),
    Knob("M4a_GMID", (5, 25)),
    Knob("M5a_GMID", (5, 25)),
    Knob("MTP_GMID", (5, 25)),
    Knob("MTN_GMID", (5, 25)),
    Knob("MP1a_L", (LMIN, LMAX)),
    Knob("MN1a_L", (LMIN, LMAX)),
    Knob("M2a_L", (LMIN, LMAX)),
    Knob("M3a_L", (LMIN, LMAX)),
    Knob("M4a_L", (LMIN, LMAX)),
    Knob("M5a_L", (LMIN, LMAX)),
    Knob("MTP_L", (LMIN, LMAX)),
    Knob("MTN_L", (LMIN, LMAX)),
    Knob("M5a_ID_over_MP1a_ID", (0.5, 4)),
    Knob("MN1a_ID_over_MP1a_ID", (0.5, 2)),
    Knob("MP1a_ID", (5e-06, 50e-06)),
    Knob("M5a_VDSAT_MARGIN", (0.0, 0.2), recorner_bound=(-0.5, 0.5)),
    Knob("M2a_VDSAT_MARGIN", (0.0, 0.2), recorner_bound=(-0.5, 0.5)),
]

# EDIT: targets. Keys must be returned by Circuit.specs().
TARGET_SPECS = {
    "GBW": Spec(5e6, "max", 1.0),
    "AC Gain (dB)": Spec(50.0, "max", 1.0),
    "PM": Spec(70.0, "max", 1.0),
    "DC CMR (dB)": Spec(50.0, "max", 0.5),
    "Area": Spec(50e-12, "min", 0.3),
    "Itotal": Spec(80e-6, "min", 1.0),
    "VOUT_DC": Spec(0.6, "eq", 2.0),
    "Output_Swing": Spec(0.48, "max", 1.0),


    "VIN_MIN": Spec(0.10, "min", 1.0),
    "VIN_MAX": Spec(1.10, "max", 1.0),
    "ICMR_overlap": Spec(0.05, "max", 1.0),
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
        maxiter=200,
        seed=1,
        n_restarts=1,
    )
