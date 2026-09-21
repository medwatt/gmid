from design import Circuit, run
from mosplot.optimizer import Corner, Knob, Spec

LUT_DIR = '/home/medwatt/coding/gmid_lookup'
NMOS, PMOS = 'nch_lvt', 'pch_lvt'

# EDIT: nominal operating conditions for this PDK.
COND = dict(vdd=1.2, vin_cm=0.6, vout_dc=0.6, cout=5e-12)

PARAMETERS = [
    Knob('M1_GMID', (5, 20)),
    Knob('M2_GMID', (10, 20)),
    Knob('M1_L', (1e-07, 1e-05)),
    Knob('M2_L', (1e-07, 1e-05)),
    Knob('M1_ID', (5e-06, 3e-05)),
    # output bias: pinned at tt; re-solved on other corners, inside the supply
    Knob("VOUT_Q", (COND["vout_dc"], COND["vout_dc"]), recorner_bound=(0.01, COND["vdd"] - 0.01)),
]

# EDIT: targets. Keys must be returned by Circuit.specs().
TARGET_SPECS = {
    "GBW": Spec(5e6, "max", 2.0),
    "AC Gain (dB)": Spec(20.0, "max", 1.0),
    "PM": Spec(70.0, "max", 1.0),
    "Area": Spec(50e-12, "min", 1.0),
    "Itotal": Spec(30e-6, "min", 1.0),
    "Output_Swing": Spec(0.6 * COND["vdd"], "max", 1.0),
    # Active-loaded output is a high-Z node; this computed VOUT_DC is what the
    # testbench sees when it applies vin_cm.
    "VOUT_DC": Spec(COND["vout_dc"], "eq", 2.0),
}

# EDIT: one Corner per process LUT (voltage corners reuse a LUT with overridden COND).
CORNERS = [
    Corner("tt", f"{LUT_DIR}/tsmc_65_lv_spectre_tt.npz", NMOS, PMOS, conditions=COND),
    Corner("ff", f"{LUT_DIR}/tsmc_65_lv_spectre_ff.npz", NMOS, PMOS, conditions=COND),
    Corner("ss", f"{LUT_DIR}/tsmc_65_lv_spectre_ss.npz", NMOS, PMOS, conditions=COND),
]

if __name__ == "__main__":
    run(workers=4, circuit=Circuit, parameters=PARAMETERS, target_specs=TARGET_SPECS,
        corners=CORNERS, output_module="./design.scs",
        maxiter=250, seed=1, n_restarts=1)
