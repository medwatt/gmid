from design import Circuit, run
from mosplot.optimizer import Corner, Knob, Spec

# EDIT: PDK lookup tables and device names.
LUT_DIR = "/home/medwatt/coding/gmid_lookup"
NMOS, PMOS = "nch_lvt", "pch_lvt"
LMIN, LMAX = 200e-9, 7.5e-6

# EDIT: nominal operating conditions for this PDK.
COND = dict(vdd=1.2, vout_dc=0.6, cout=1e-12, iref=5e-6, k=4.0)

# EDIT: knob bounds. Names must match the circuit's KNOBS.
PARAMETERS = [
    Knob("M1_GMID", (8.0, 15.0)),
    Knob("M1_L", (LMIN, LMAX)),
    Knob("M3_GMID", (8.0, 15.0)),
    Knob("M3_L", (LMIN, LMAX)),
    Knob("VDSAT_MARGIN", (0.03, 0.2), recorner_bound=(-0.5, 0.5)),
]

# EDIT: targets. Keys must be returned by Circuit.specs().
TARGET_SPECS = {
    "Rout": Spec(5e9, "max", 1.0),
    "Area": Spec(20e-12, "min", 10.0),
    "Vcompliance": Spec(0.2, "min", 2.0),
    "Margin min":   Spec(0.02, "max", 1.0, scale=0.05),
}

# EDIT: one Corner per process LUT (voltage corners reuse a LUT with overridden COND).
CORNERS = [
    Corner("tt", f"{LUT_DIR}/tsmc_65_lv_spectre_tt.npz", NMOS, PMOS, conditions=COND),
    Corner("ss", f"{LUT_DIR}/tsmc_65_lv_spectre_ss.npz", NMOS, PMOS, conditions=COND),
    Corner("ff", f"{LUT_DIR}/tsmc_65_lv_spectre_ff.npz", NMOS, PMOS, conditions=COND),
]

if __name__ == "__main__":
    run(
        workers=4,
        circuit=Circuit,
        parameters=PARAMETERS,
        target_specs=TARGET_SPECS,
        corners=CORNERS,
        output_module="./wide_swing_pmos_mirror.scs",
        maxiter=500,
        seed=1,
        n_restarts=1,
    )
