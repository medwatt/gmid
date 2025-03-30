# MOSFET Characterization in Python

## Introduction

This tool does three things:

1. Generate lookup tables of MOSFET parameters.

2. Make plots of MOSFET parameters.

3. Optimize analog circuits.


## Installation

### Requirements

- Python (3.9+ recommended)

- Numpy, Scipy, Matplotlib

Also install a SPICE simulator such as [ngspice](https://ngspice.sourceforge.io/) or HSPICE to generate lookup tables.

### Steps

- Clone the repository:

```bash
git clone https://github.com/medwatt/gmid.git
```

- Install the package:

```bash
pip install .
```

## Generating a Lookup Table

Before plotting, create a lookup table with all the MOSFET parameters. For example:

```python
from mosplot.lookup_table_generator import LookupTableGenerator

obj = LookupTableGenerator(
    description="freepdk 45nm",

    # Simulator to use
    simulator="ngspice", # "ngspice" or "hspice"

    # Provide path to simulator if not in system path
    simulator_path="/usr/bin/ngspice",

    # Files to include with `.INCLUDE`
    include_paths=[
        "/home/username/gmid/models/NMOS_VTH.lib",
        "/home/username/gmid/models/PMOS_VTH.lib",
    ],

    # Names of models to simulate
    model_names={
        "NMOS_VTH": "nmos",
        "PMOS_VTH": "pmos",
    },

    # Symbols for detecting the transistor:
    mos_spice_symbols = ("m1", "m1"),

    # Fixed width for all simulations
    width=10e-6,

    # Voltage sweep parameters
    vsb=(0, 1.0, 0.1),
    vgs=(0, 1.0, 0.01),
    vds=(0, 1.0, 0.01),

    # Length sweep parameters
    lengths=[50e-9, 100e-9, 200e-9, 400e-9, 800e-9, 1.6e-6, 3.2e-6, 6.4e-6],
)

# Optionally, run an op simulation to check outputs.
# obj.op_simulation()

# Build and store the table
obj.build("./freepdk_45nm")
```

## Using the Tool

It is recommended to use this tool in a Jupyter notebook.

### Imports

```python
import numpy as np
from mosplot.plot import load_lookup_table, Mosfet
```

Load a lookup table:

```python
lookup_table = load_lookup_table("path/to/lookup-table.npy")
```

The `Mosfet` class helps you make plots. If you want to change the style,
import `matplotlib`:

```python
import matplotlib.pyplot as plt
plt.style.use('path/to/style')
```

## Making Simple Plots

Create a MOSFET instance:

```python
nmos = Mosfet(lookup_table=lookup_table, mos="NMOS_VTH", vsb=0.0, vds=0.5, vgs=(0.3, 1))
```

The above code filters the table at `vsb=0.0` and `vds=0.5` for all lengths
and for `vgs` values between 0.3 and 1 (or add a step like `(0.3, 1,
0.02)`). To select all values of `vgs`, set it to `None` or omit it.

To create a plot, use `plot_by_expression`. For example, to plot $I_{D}/W$
vs $g_{m}/I_{D}$:

```python
nmos.plot_by_expression(
    x_expression = nmos.gmid_expression,
    y_expression = nmos.current_density_expression,
)
```

![current density plot](./figures/nmos_current_density.svg)

Some common expressions include:

- `gmid_expression`

- `vgs_expression`

- `vds_expression`

- `vsb_expression`

- `gain_expression`

- `current_density_expression`

- `transist_frequency_expression`

- `early_voltage_expression`

You can pass a list of lengths to the `lengths` parameter:

```python
nmos.plot_by_expression(
    x_expression = nmos.gmid_expression,
    y_expression = nmos.current_density_expression,
    lengths = [5.0e-08, 1.0e-07, 2.0e-07]
)
```

![current density plot](./figures/nmos_current_density_filtered.svg)

The tool tries to choose proper axis scales. For example, the last plot may use
a `log` scale. You can override options:

```python
nmos.plot_by_expression(
    x_expression = nmos.gmid_expression,
    y_expression = nmos.current_density_expression,
    lengths = [5.0e-08, 1.0e-07, 2.0e-07],
    y_scale = 'linear',
    x_limit = (5, 20),
    y_limit = (0, 300),
    save_fig="path/to/save/figure/with/extension"
)
```

![current density plot](./figures/nmos_current_density_options.svg)

To plot something custom, for example $V_{GS}$ on the x-axis and a custom
expression on the y-axis:

```python
from mosplot.plot import Expression

nmos.plot_by_expression(
    x_expression = nmos.vgs_expression,
    y_expression = Expression(
        variables=["id", "gds"],
        function=lambda x, y: x / y,
        label="$I_D / g_{ds} (A/S)$",
    )
)
```

![custom expression](./figures/nmos_custom_expression_1.svg)

## Looking Up Values

You can also get raw values by interpolation. For example, to look up the gain
at a specific point:
```python
x = nmos.interpolate(
    x_expression=nmos.lengths_expression,
    x_value=100e-9,
    y_expression=nmos.gmid_expression,
    y_value=15,
    z_expression=nmos.gain_expression,
)
```

Or to get a 2D sweep:

```python
x = nmos.interpolate(
    x_expression=nmos.vdsat_expression,
    x_value=(0.08, 0.12, 0.01),
    y_expression=nmos.gds_expression,
    y_value=(1e-6, 4e-6, 1e-6),
    z_expression=nmos.gmid_expression,
)
```

To look up an expression without interpolation:

```python
x = nmos.lookup_expression_from_table(
    lengths=100e-9,
    vsb=0,
    vds=(0.0, 1, 0.01),
    vgs=(0.0, 1.01, 0.2),
    primary="vds",
    expression=nmos.current_density_expression,
)
```

## Plotting Methods

### Plot by Sweep

Use `plot_by_sweep` to create various plots. For example, to plot the output
characteristic of a MOSFET:

```python
nmos.plot_by_sweep(
    lengths=180e-9,
    vsb = 0,
    vds = (0.0, 1, 0.01), # you can also set to `None`
    vgs = (0.0, 1.01, 0.2),
    x_expression = nmos.vds_expression,
    y_expression = nmos.id_expression,
    primary = "vds",
    x_eng_format=True,
    y_eng_format=True,
    y_scale='linear',
)
```

![output characteristic](./figures/nmos_output_characteristics.svg)

### Quick Plot

Combine data from multiple plots into one. For example:

```python
vdsat = nmos.plot_by_expression(
    lengths=[45e-9],
    x_expression = nmos.vgs_expression,
    y_expression = nmos.vdsat_expression,
    return_result = True,
)

vov = nmos.plot_by_expression(
    lengths=[45e-9],
    x_expression = nmos.vgs_expression,
    y_expression = Expression(
        variables=["vgs", "vth"],
        function=lambda x, y: x - y,
    ),
    return_result = True,
)

vstar = nmos.plot_by_expression(
    lengths=[45e-9],
    x_expression = nmos.vgs_expression,
    y_expression = Expression(
        variables=["gm", "id"],
        function=lambda x, y: 2 / (x/y),
    ),
    return_result=True,
)
```

Then, combine the results with `quick_plot`:

```python
nmos.quick_plot(
    x = [vdsat[0][0], vstar[0][0], vov[0][0]],
    y = [vdsat[1][0], vstar[1][0], vov[1][0]],
    legend = ["$V_{\\mathrm{DS}_{\\mathrm{SAT}}}$", "$V^{\\star}$", "$V_{\\mathrm{OV}}$"],
    x_limit = (0.1, 1),
    y_limit = (0, 0.6),
    x_label = "$V_{\\mathrm{GS}}$",
    y_label = "$V$",
)
```

![qucik plot](./figures/nmos_quick_plot.svg)

## Optimization

You can optimize analog circuits using this tool. Define optimization
parameters and target specifications, then run the optimizer.

```python
from mosplot.optimizer import Optimizer, DesignReport
from datatypes import Spec, OptimizationParameter

parameters = [
    OptimizationParameter("L_input", (100e-9, 2e-6)),
    OptimizationParameter("gmid_input", (7, 16)),
    OptimizationParameter("L_load", (100e-9, 2e-6)),
    OptimizationParameter("gmid_load", (7, 15)),
    OptimizationParameter("L_tail", (100e-9, 2e-6)),
    OptimizationParameter("gmid_tail", (7, 15)),
    OptimizationParameter("Ibias", (1e-6, 40e-6))
]

target_specs = {
    "GBW":       Spec(5e6, "max", 5),
    "Gain":      Spec(30, "max", 2),
    "Ibias":     Spec(20e-6, "min", 1),
    "ICMR_LOW":  Spec(0.1, "min", 1),
    "ICMR_HIGH": Spec(0.7, "max", 5),
    "CMRR":      Spec(2000, "max", 3),
    "Area":      Spec(15e-12, "min", 1)
}
```

Create your circuit instance that defines how the specs are computed and run
the optimizer:

```python
circuit = Circuit(lookup_table, pmos_range, nmos_range, "pch_lvt", "nch_lvt", VDD, CL)
optimizer = Optimizer(circuit, parameters, target_specs)
optimizer.optimize(maxiter=5)

report = DesignReport(circuit, optimizer)
print("\nDesign Report:")
print(report.report())
```

## Acknowledgment
- HSPICE output parsing is based on [this script](https://github.com/HMC-ACE/hspiceParser) .

- If you find this tool useful, please cite it.
