from design import Circuit, run
from mosplot.optimizer import Corner, Knob, Spec

LUT_DIR = "/home/medwatt/coding/gmid_lookup"
NMOS, PMOS = "nch_lvt", "pch_lvt"
COND = dict(vdd=1.2)

PARAMETERS = [
    Knob("M2_GMID", (10.0, 20.0)),
    Knob("M4_GMID", (8.0, 18.0)),
    Knob("M5_GMID", (8.0, 18.0)),
    Knob("M7_GMID", (10.0, 20.0)),
    Knob("L_n", (500e-9, 4e-6)),
    Knob("L_nc", (300e-9, 4e-6)),
    Knob("L_pc", (300e-9, 4e-6)),
    Knob("L_pm", (500e-9, 4e-6)),
    Knob("K", (4.0, 8.0)),
    Knob("IREF", (10e-6, 80e-6)),
]

TARGET_SPECS = {
    "gm_ref": Spec(150e-6, "max", 5.0),
    "gmR": Spec(1.0, "eq", 20.0),
    "Area": Spec(80e-12, "min", 0.5),
    "Itotal": Spec(60e-6, "min", 0.5),
}

CORNERS = [
    Corner("tt", f"{LUT_DIR}/tsmc_65_lv_spectre_tt.npz", NMOS, PMOS, conditions=COND),
    Corner("ff", f"{LUT_DIR}/tsmc_65_lv_spectre_ff.npz", NMOS, PMOS, conditions=COND),
    Corner("ss", f"{LUT_DIR}/tsmc_65_lv_spectre_ss.npz", NMOS, PMOS, conditions=COND),
    Corner("tt_14", f"{LUT_DIR}/tsmc_65_lv_spectre_tt.npz", NMOS, PMOS, conditions=dict(vdd=1.4)),
    Corner("ff_14", f"{LUT_DIR}/tsmc_65_lv_spectre_ff.npz", NMOS, PMOS, conditions=dict(vdd=1.4)),
    Corner("ss_14", f"{LUT_DIR}/tsmc_65_lv_spectre_ss.npz", NMOS, PMOS, conditions=dict(vdd=1.4)),
    Corner("tt_10", f"{LUT_DIR}/tsmc_65_lv_spectre_tt.npz", NMOS, PMOS, conditions=dict(vdd=1.0)),
    Corner("ff_10", f"{LUT_DIR}/tsmc_65_lv_spectre_ff.npz", NMOS, PMOS, conditions=dict(vdd=1.0)),
    Corner("ss_10", f"{LUT_DIR}/tsmc_65_lv_spectre_ss.npz", NMOS, PMOS, conditions=dict(vdd=1.0)),
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
