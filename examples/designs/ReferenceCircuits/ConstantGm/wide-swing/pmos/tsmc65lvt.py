from design import Circuit, run
from mosplot.optimizer import Corner, Knob, Spec

LUT_DIR = "/home/medwatt/coding/gmid_lookup"
NMOS, PMOS = "nch_lvt", "pch_lvt"
COND = dict(vdd=1.2)

PARAMETERS = [
    Knob("M2_GMID", (10.0, 20.0)),
    Knob("M3_GMID", (8.0, 18.0)),
    Knob("M5_GMID", (8.0, 18.0)),
    Knob("M7_GMID", (10.0, 20.0)),
    Knob("L_p",  (500e-9, 7.5e-6)),
    Knob("L_pc", (500e-9, 7.5e-6)),
    Knob("L_nc", (500e-9, 7.5e-6)),
    Knob("L_n",  (500e-9, 7.5e-6)),
    Knob("K",    (4.0, 8.0)),
    Knob("IREF", (10e-6, 100e-6)),
    Knob("Rn_V", (0.05, 0.5)),
    Knob("Rp_V", (0.05, 0.5)),
]

TARGET_SPECS = {
    "gm_ref":     Spec(150e-6, "max", 5.0),
    "gmR":        Spec(1.0, "eq", 20.0),
    "VDS_margin": Spec(0.10, "max", 10.0),
    "Area":       Spec(150e-12, "min", 0.5),
    "Itotal":     Spec(80e-6,  "min", 0.5),
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
        output_module="./design.scs",
        maxiter=80,
        seed=1,
        n_restarts=1,
    )
