from design import Circuit, run
from mosplot.optimizer import Corner, Knob, Spec

LUT_DIR = "/home/medwatt/coding/gmid_lookup"
NMOS, PMOS = "nch_lvt", "pch_lvt"
LMIN, LMAX = 200e-9, 7.5e-6

COND = dict(vdd=1.2, vout_dc=0.6, cout=1e-12, iref=5e-6, k=4.0)

PARAMETERS = [
    Knob("M1_GMID", (8.0, 15.0)),
    Knob("M1_L", (LMIN, LMAX)),
    Knob("M3_GMID", (8.0, 15.0)),
    Knob("M3_L", (LMIN, LMAX)),
]

TARGET_SPECS = {
    "Rout": Spec(100e6, "max", 2.0),
    "Area": Spec(20e-12, "min", 10.0),
    "Vcompliance": Spec(0.3, "min", 2.0),
}

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
        output_module="./wilson_nmos_mirror.scs",
        maxiter=500,
        seed=1,
        n_restarts=1,
    )
