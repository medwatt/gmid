import numpy as np
from mosplot.plot import Mosfet, Expression, load_lookup_table
from mosplot.optimizer import Optimizer, DesignReport, Spec, OptimizationParameter


class MillerOpampCircuit:
    def __init__(self, lookup_table_path, pmos_range, nmos_range, pmos_name, nmos_name, VDD, CL):
        self.lookup_table = load_lookup_table(lookup_table_path)
        self.pmos_range = pmos_range
        self.nmos_range = nmos_range
        self.pmos_name = pmos_name
        self.nmos_name = nmos_name
        self.VDD = VDD
        self.CL = CL
        self.gdsid_expression = Expression(variables=["gds", "id"], function=lambda gds, id: gds / id)
        self.map_transistors()

    def map_transistors(self):
        # first stage: 3‐device stack => |vds|=VDD/3
        self.pmos_s1 = Mosfet(
            lookup_table=self.lookup_table,
            mos=self.pmos_name,
            vbs=0.0,
            vds=-self.VDD / 3,
            vgs=self.pmos_range,
        )
        self.nmos_s1 = Mosfet(
            lookup_table=self.lookup_table,
            mos=self.nmos_name,
            vbs=0.0,
            vds=self.VDD / 3,
            vgs=self.nmos_range,
        )
        # second stage: 2‐device stack => |vds|=VDD/2
        self.pmos_s2 = Mosfet(
            lookup_table=self.lookup_table,
            mos=self.pmos_name,
            vbs=0.0,
            vds=-self.VDD / 2,
            vgs=self.pmos_range,
        )
        self.nmos_s2 = Mosfet(
            lookup_table=self.lookup_table,
            mos=self.nmos_name,
            vbs=0.0,
            vds=self.VDD / 2,
            vgs=self.nmos_range,
        )

        self.transistor_mapping = {
            ########################################
            #         Stage 1 Transistors          #
            ########################################
            "pmos_s1_input": {
                "transistor": self.pmos_s1,
                "length": "L_s1_input",
                "gmid": "gmid_s1_input",
                "current": None,
            },
            "nmos_s1_load": {
                "transistor": self.nmos_s1,
                "length": "L_s1_load",
                "gmid": "gmid_s1_load",
                "current": None,
            },
            "pmos_s1_tail": {
                "transistor": self.pmos_s1,
                "length": "L_s1_tail",
                "gmid": "gmid_s1_tail",
                "current": None,
            },
            ########################################
            #         Stage 2 Transistors          #
            ########################################
            # Load transistor of the second stage is part of a current-mirror configuration.
            # Its length and vsg are determined by the length and vsg of first stage's tail transistor.
            "pmos_s2_load": {
                "transistor": self.pmos_s2,
                "length": "L_s1_tail",
                "gmid": "gmid_s2_load",
                "current": None,
            },
            "nmos_s2_input": {
                "transistor": self.nmos_s2,
                "length": "L_s2_input",
                "gmid": "gmid_s2_input",
                "current": None,
            },
        }

    def compute_device_dimensions(self, params):
        dims = {}
        for key, mapping in self.transistor_mapping.items():
            transistor = mapping["transistor"]
            length = params[mapping["length"]]
            gmid = params[mapping["gmid"]]
            current = mapping["current"]
            idw = transistor.interpolate(
                x_expression=transistor.length_expression,
                x_value=length,
                y_expression=transistor.gmid_expression,
                y_value=gmid,
                z_expression=transistor.current_density_expression,
                fast=True,
            )
            width = (current / idw)[0]
            area = width * length
            dims[key] = {
                "Length": length,
                "Width": width,
                "Area": area,
                "Current": current,
                "GMID": gmid,
            }
        self.device_dimensions = dims

    def evaluate_specs(self, **params):
        # Assign currents.
        s1_tail_current = params["Is1"]
        s2_input_current = params["Is2"]
        s1_input_current = s1_tail_current / 2.0
        self.transistor_mapping["pmos_s1_input"]["current"] = s1_input_current
        self.transistor_mapping["nmos_s1_load"]["current"] = s1_input_current
        self.transistor_mapping["pmos_s1_tail"]["current"] = s1_tail_current
        self.transistor_mapping["nmos_s2_input"]["current"] = s2_input_current
        self.transistor_mapping["pmos_s2_load"]["current"] = s2_input_current

        # Compensation capacitor.
        comp_cap = params["Cc"]

        ################################################################################
        #                              Stage 1 Parameters                              #
        ################################################################################
        cdd_s1_input, vsg_s1_input, vdsat_s1_input, gdsid_s1_input = self.pmos_s1.interpolate(
            x_expression=self.pmos_s1.length_expression,
            x_value=params["L_s1_input"],
            y_expression=self.pmos_s1.gmid_expression,
            y_value=params["gmid_s1_input"],
            z_expression=[
                self.pmos_s1.cdd_expression,
                self.pmos_s1.vsg_expression,
                self.pmos_s1.vdsat_expression,
                self.gdsid_expression,
            ],
            fast=True,
        )
        cdd_s1_load, vgs_s1_load, gdsid_s1_load = self.nmos_s1.interpolate(
            x_expression=self.nmos_s1.length_expression,
            x_value=params["L_s1_load"],
            y_expression=self.nmos_s1.gmid_expression,
            y_value=params["gmid_s1_load"],
            z_expression=[
                self.nmos_s1.cdd_expression,
                self.nmos_s1.vgs_expression,
                self.gdsid_expression,
            ],
            fast=True,
        )
        vgs_s1_tail, vdsat_s1_tail, early_s1_tail = self.pmos_s1.interpolate(
            x_expression=self.pmos_s1.length_expression,
            x_value=params["L_s1_tail"],
            y_expression=self.pmos_s1.gmid_expression,
            y_value=params["gmid_s1_tail"],
            z_expression=[
                self.pmos_s1.vgs_expression,
                self.pmos_s1.vdsat_expression,
                self.pmos_s1.early_voltage_expression,
            ],
            fast=True,
        )

        ################################################################################
        #                              Stage 2 Parameters                              #
        ################################################################################
        # The vsg of the load transistor of stage 2 is determined by the vsg of
        # the tail transistor of the first stage.
        gdsid_s2_load, gmid_s2_load = self.pmos_s2.interpolate(
            x_expression=self.pmos_s2.length_expression,
            x_value=params["L_s1_tail"],
            y_expression=self.pmos_s2.vgs_expression,
            y_value=vgs_s1_tail[0],
            z_expression=[
                self.gdsid_expression,
                self.pmos_s2.gmid_expression,
            ],
            fast=True,
        )
        # Determine gmid from length and vsg
        params["gmid_s2_load"] = gmid_s2_load[0]

        cgs_s2_input, gdsid_s2_input = self.nmos_s2.interpolate(
            x_expression=self.nmos_s2.length_expression,
            x_value=params["L_s2_input"],
            y_expression=self.nmos_s2.gmid_expression,
            y_value=params["gmid_s2_input"],
            z_expression=[
                self.nmos_s2.cgs_expression,
                self.gdsid_expression
            ],
            fast=True,
        )

        ################################################################################
        #                              Specs Computation                               #
        ################################################################################
        # GBW
        gm1 = params["gmid_s1_input"] * s1_input_current
        gbw = gm1 / (2 * np.pi * (comp_cap + cdd_s1_input + cdd_s1_load + cgs_s2_input)[0])

        # DC gain
        gain_s1 = params["gmid_s1_input"] / (gdsid_s1_input[0] + gdsid_s1_load[0])
        gain_s2 = params["gmid_s2_input"] / (gdsid_s2_input[0] + gdsid_s2_load[0])
        gain = gain_s1 * gain_s2

        # ICMR low, high
        icmr_low = vgs_s1_load[0] - vsg_s1_input[0] + vdsat_s1_input[0]
        icmr_high = self.VDD - vsg_s1_input[0] - vdsat_s1_tail[0]

        # CMRR
        cmrr = gain_s1 * params["gmid_s1_load"] * early_s1_tail[0]

        # Slew Rate
        sr = s1_tail_current / (comp_cap + cdd_s1_input + cdd_s1_load + cgs_s2_input)[0]

        # Phase Margin
        k = params["gmid_s1_input"] / params["gmid_s2_input"]
        a1 = k * (s1_input_current / s2_input_current)
        a2 = a1 * (self.CL / comp_cap)
        pm = 90 - np.degrees(np.arctan(a1)) - np.degrees(np.arctan(a2))

        # Total Current
        itotal = s1_tail_current + s2_input_current

        # Area
        self.compute_device_dimensions(params)
        area = sum(d["Area"] for d in self.device_dimensions.values())

        return {
            "GBW": gbw,
            "Gain": gain,
            "Ibias": s1_tail_current,
            "Itotal": itotal,
            "ICMR_LOW": icmr_low,
            "ICMR_HIGH": icmr_high,
            "CMRR": cmrr,
            "PhaseMargin": pm,
            "SlewRate": sr,
            "Area": area,
        }


if __name__ == "__main__":
    lookup = "/home/medwatt/git/gmid/tests/lookup_table_generator/test/tsmc65_test/tsmc65lvt.npz"
    VDD = 1.2
    CL = 5e-12
    pmos_range = (-VDD, -0.1)
    nmos_range = (0.1, VDD)
    pmos_name, nmos_name = "pch_lvt", "nch_lvt"

    circuit = MillerOpampCircuit(
        lookup, pmos_range, nmos_range, pmos_name, nmos_name, VDD, CL
    )

    parameters = [
        OptimizationParameter("L_s1_input", (100e-9, 3e-6)),
        OptimizationParameter("gmid_s1_input", (7, 20)),
        OptimizationParameter("L_s1_load", (100e-9, 3e-6)),
        OptimizationParameter("gmid_s1_load", (7, 15)),
        OptimizationParameter("L_s1_tail", (100e-9, 3e-6)),
        OptimizationParameter("gmid_s1_tail", (7, 20)),
        OptimizationParameter("L_s2_input", (100e-9, 3e-6)),
        OptimizationParameter("gmid_s2_input", (7, 20)),
        OptimizationParameter("Is1", (10e-6, 40e-6)),
        OptimizationParameter("Is2", (20e-6, 100e-6)),
        OptimizationParameter("Cc", (0.1e-12, 5e-12)),
    ]

    target_specs = {
        "GBW": Spec(5e6, "max", 3),
        "Itotal": Spec(80e-6, "min", 3),
        "Gain": Spec(10 ** (60 / 20), "max", 3),
        "ICMR_LOW": Spec(0.1, "min", 1),
        "ICMR_HIGH": Spec(0.6, "max", 1),
        "CMRR": Spec(10 ** (70 / 20), "max", 1),
        "PhaseMargin": Spec(70, "max", 5),
        "SlewRate": Spec(15e6, "max", 2),
        "Area": Spec(20e-12, "min", 0.1),
    }

    opt = Optimizer(circuit, parameters, target_specs)
    opt.optimize(maxiter=30)

    report = DesignReport(circuit, opt)
    print(report.report())

    best = opt.get_opt_params()

    vsg1 = circuit.pmos_s1.interpolate(
        x_expression=circuit.pmos_s1.length_expression,
        x_value=best["L_s1_input"],
        y_expression=circuit.pmos_s1.gmid_expression,
        y_value=best["gmid_s1_input"],
        z_expression=circuit.pmos_s1.vsg_expression,
    )

    print("")
    print(f"pmos_s1_input VSG = {vsg1[0]}")
    print(f"Cc = {best["Cc"]}")
