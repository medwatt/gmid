import numpy as np
from mosplot.plot import Mosfet, Expression, load_lookup_table
from mosplot.optimizer import Optimizer, DesignReport, Spec, OptimizationParameter


# Circuit <<<
class Circuit:
    def __init__(self, lookup_table_path, pmos_range, nmos_range, pmos_name, nmos_name, VDD, CL):
        self.lookup_table = load_lookup_table(lookup_table_path)
        self.pmos_range = pmos_range
        self.nmos_range = nmos_range
        self.pmos_name = pmos_name
        self.nmos_name = nmos_name
        self.VDD = VDD
        self.CL = CL
        self.gdsid_expression = Expression(variables=["gds", "id"], function=lambda x, y: x / y)
        self.map_transistors()

    def map_transistors(self):
        self.pmos = Mosfet(
            lookup_table=self.lookup_table,
            mos=self.pmos_name,
            vsb=0.0,
            vds=-self.VDD / 3,
            vgs=self.pmos_range
        )
        self.nmos = Mosfet(
            lookup_table=self.lookup_table,
            mos=self.nmos_name,
            vsb=0.0,
            vds=self.VDD / 3,
            vgs=self.nmos_range,
        )

        # We have 3 kinds of transistors.
        self.transistor_mapping = {
            "pmos_input": {
                "transistor": self.pmos,
                "length": "L_input",
                "gmid": "gmid_input",
                "current": None,
            },
            "nmos_load": {
                "transistor": self.nmos,
                "length": "L_load",
                "gmid": "gmid_load",
                "current": None,
            },
            "pmos_tail": {
                "transistor": self.pmos,
                "length": "L_tail",
                "gmid": "gmid_tail",
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
        L_input    = params["L_input"]
        gmid_input = params["gmid_input"]
        L_load     = params["L_load"]
        gmid_load  = params["gmid_load"]
        L_tail     = params["L_tail"]
        gmid_tail  = params["gmid_tail"]
        Ibias      = params["Ibias"]

        self.transistor_mapping["pmos_input"]["current"] = Ibias / 2.0
        self.transistor_mapping["nmos_load"]["current"] = Ibias / 2.0
        self.transistor_mapping["pmos_tail"]["current"] = Ibias

        cdd_input, vgs_input, vdsat_input, gdsid_input = self.pmos.interpolate(
            x_expression=self.pmos.length_expression,
            x_value=L_input,
            y_expression=self.pmos.gmid_expression,
            y_value=gmid_input,
            z_expression=[
                self.pmos.cdd_expression,
                self.pmos.vgs_expression,
                self.pmos.vdsat_expression,
                self.gdsid_expression,
            ],
            fast=True,
        )
        cdd_load, vgs_load, gdsid_load = self.nmos.interpolate(
            x_expression=self.nmos.length_expression,
            x_value=L_load,
            y_expression=self.nmos.gmid_expression,
            y_value=gmid_load,
            z_expression=[
                self.nmos.cdd_expression,
                self.nmos.vgs_expression,
                self.gdsid_expression,
            ],
            fast=True,
        )
        vdsat_tail, gdsid_tail = self.pmos.interpolate(
            x_expression=self.pmos.length_expression,
            x_value=L_tail,
            y_expression=self.pmos.gmid_expression,
            y_value=gmid_tail,
            z_expression=[
                self.pmos.vdsat_expression,
                self.gdsid_expression,
            ],
            fast=True,
        )

        current_input = Ibias / 2.0
        gbw = (gmid_input * current_input) / (2 * np.pi * (self.CL + cdd_input[0] + cdd_load[0]))
        gain = (gmid_input / (gdsid_input + gdsid_load))[0]
        icmr_low = (vgs_load + vgs_input + vdsat_input)[0]
        icmr_high = (self.VDD + vgs_input - vdsat_tail)[0]
        cmrr = (gain * gmid_load / gdsid_tail)[0]

        self.compute_device_dimensions(params)
        overall_area = sum(item["Area"] for item in self.device_dimensions.values() if item["Area"])

        specs = {
            "GBW": gbw,
            "Gain": gain,
            "Ibias": Ibias,
            "ICMR_LOW": icmr_low,
            "ICMR_HIGH": icmr_high,
            "CMRR": cmrr,
            "Area": overall_area,
        }

        return specs
# >>>

if __name__ == "__main__":
    # Define circuit constants and lookup table location.
    lookup_table = "/home/medwatt/git/gmid/tests/lookup_table_generator/tsmc65/tsmc_65nm.npy"
    VDD = 1.2
    CL = 5e-12
    pmos_range = (-1.2, -0.1)
    nmos_range = (0.1, 1.2)
    pmos_name = "pch_lvt"
    nmos_name = "nch_lvt"

    # Create the circuit instance.
    circuit = Circuit(lookup_table, pmos_range, nmos_range, pmos_name, nmos_name, VDD, CL)

    # Define optimization parameters as a list of OptimizationParameter objects.
    parameters = [
        OptimizationParameter("L_input", (100e-9, 2e-6)),
        OptimizationParameter("gmid_input", (7, 16)),
        OptimizationParameter("L_load", (100e-9, 2e-6)),
        OptimizationParameter("gmid_load", (7, 15)),
        OptimizationParameter("L_tail", (100e-9, 2e-6)),
        OptimizationParameter("gmid_tail", (7, 15)),
        OptimizationParameter("Ibias", (1e-6, 40e-6))
    ]

    # Define target specifications using Spec objects.
    target_specs = {
        "GBW":       Spec(5e6, "max", 5),
        "Gain":      Spec(30, "max", 2),
        "Ibias":     Spec(20e-6, "min", 1),
        "ICMR_LOW":  Spec(0.1, "min", 1),
        "ICMR_HIGH": Spec(0.7, "max", 5),
        "CMRR":      Spec(2000, "max", 3),
        "Area":      Spec(15e-12, "min", 1)
    }

    # Instantiate and run the optimizer.
    optimizer = Optimizer(circuit, parameters, target_specs)
    result = optimizer.optimize(maxiter=5)

    # Generate and print the design report.
    report = DesignReport(circuit, optimizer)
    print("\nDesign Report:")
    print(report.report())

    # Get optimal parameters for interpolation.
    opt_params = optimizer.get_opt_params()

    # Interpolate VGS for the input device and tail.
    vgs_input = circuit.pmos.interpolate(
        x_expression=circuit.pmos.length_expression,
        x_value=opt_params["L_input"],
        y_expression=circuit.pmos.gmid_expression,
        y_value=opt_params["gmid_input"],
        z_expression=circuit.pmos.vgs_expression
    )
    vgs_tail = circuit.pmos.interpolate(
        x_expression=circuit.pmos.length_expression,
        x_value=opt_params["L_tail"],
        y_expression=circuit.pmos.gmid_expression,
        y_value=opt_params["gmid_tail"],
        z_expression=circuit.pmos.vgs_expression
    )

    print("")
    print(f"Input Pair VSG = {-vgs_input[0]}")
    print(f"Tail VG = {VDD + vgs_tail[0]}")
