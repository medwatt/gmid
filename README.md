# MOSFET Characterization in Python

## Motivation

This tool has the following goals:

1. Provide an easy way of creating plots of MOSFET parameters, such as those
   used in the gm/ID design methodology.

2. Provide a tool that does not depend on any proprietary software or require
   licensing fees.

3. Open source so that it can be easily modified/extended by the user.

## Installation

### Requirements

This tools is written in Python and requires the following:

- `Numpy`, `Scipy`, and `Matplotlib` for data analysis and plotting.

-  [`ngspice`](https://ngspice.sourceforge.io/) or `hspice` for generating the
   lookup table.

## Generating a Lookup Table

Before any plots can be made, a lookup table of all the relevant parameters
must first be created. This is done by instantiating an object from the
`LookupTableGenerator` and then building the table with the `build` method. An
example is given below.

```python
from main.lookup_table_generator import LookupTableGenerator

obj = LookupTableGenerator(
    description="freepdk 45nm ngspice",
    simulator="ngspice",
    model_paths=[
        "/home/username/gmid/models/NMOS_VTH.lib",
        "/home/username/gmid/models/PMOS_VTH.lib",
        ],
    model_names={
        "nmos": "NMOS_VTH",
        "pmos": "PMOS_VTH",
    },
    vsb=(0, 1.0, 0.1),
    vgs=(0, 1.0, 0.01),
    vds=(0, 1.0, 0.01),
    width=10e-6,
    lengths=[50e-9, 100e-9, 200e-9, 400e-9, 800e-9, 1.6e-6, 3.2e-6, 6.4e-6],
)
obj.build("/home/username/gmid/lookup_tables/freepdk_45nm_ngspice.npy")
```

A summary of some of the parameters is given below:

- The simulator used is specified with the `simulator` parameter. At the
  moment, only `ngspice` and `hspice` are supported. If you're using windows or
  some linux distribution where `ngspice` and `hspice` are named differently,
  you will have to modify the `NGSPICE_PATH` and `HSPICE_PATH` variables inside
  the `lookup_table_generator.py` file to point to the binaries on your system.

- The lookup_table will be generated for a specific transistor model. Provide
  the location of the model files as a list using the `model_paths` parameter.
  Since it is possible to have more than one model definition inside a file,
  you need to specify the model name. This is done via the `model_names`
  parameter, where the keys are always `"nmos"` and `"pmos` and their values
  are the names of the models to be used.

- If there's a specific need to pass in some custom SPICE commands, these
  should be done via the `raw_spice` parameter (not shown in the example above).

- To generate a lookup table, the bulk, gate, and drain voltages relative to
  the source have to be swept over a range of voltages. Specify the range in
  the form `(start, stop, step)`. The smaller the step size, the bigger is the
  size of the lookup table.

- The `lengths` can be provided as a list of discrete values or a 1-dimensional
  `numpy` array.

- Only a single `width` should be provided. The assumption here is that the
  parameters of the MOSFET scale linearly with the width. Because of this
  assumption, all parameters that are width-dependent must be de-normalized
  with respect to the current or width that you're working with.

- The directory where the generated lookup table is saved is passed directly to
  the `build` method.

## Using the Tool

Because of the interactive nature of designing analog circuits, using this
script within a `jupyter` notebook is highly recommended.

### Imports

If the `gmid.py` file is located somewhere not in the same directory as the
file you're running, add its path so that the interpreter knows about it.

```python
import sys
sys.path.insert(0, "path/to/parent/directory/of/gmid.py file")
```

We begin by making the following imports:

```python
import numpy as np
from gmid import load_lookup_table, GMID
```


The `load_lookup_table` function loads a lookup table such as the one generated
in the previous section.

```python
lookup_table = load_lookup_table("path/to/lookup-table.npy")
```

The `GMID` class contains methods that can be used to generate plots
seamlessly. If you plan to modify the style of the plots or plot things
differently, you will also have to import `matplotlib`.

```python
import matplotlib.pyplot as plt
plt.style.use('path/to/style')
```

### Making Simple Plots

We start by creating an object called `nmos` that selects the NMOS
from the lookup table and sets the source-bulk and drain-source voltages to
some fixed values. Since the data is 4-dimensional, it is necessary to fix two
of the variables at a time to enable 2-dimensional plotting.

```python
nmos = GMID(lookup_table=lookup_table, mos="nmos", vsb=0.0, vds=0.5, vgs=(0.3, 1))
```

The above code filters the table at `vsb=0.0` and `vds=0.5` for all `lengths`
and for `vgs` values between `(0.3, 1)`. You can also include a step such as
`(0.3, 1, 0.02)`. If you want all values of `vgs`, either set it to `None` or
don't include it.

Methods are available to create the most commonly-used plots in the gm/ID
methodology so that you don't have to type them. These are:

- `current_density_plot()`: this plots $I_{D}/W$ vs $g_{m}/I_{D}$.
- `gain_plot()`: this plots $g_m / g_{ds}$ vs $g_{m}/I_{D}$.
- `transit_frequency_plot()`: this plots $f_{T}$ vs $g_{m}/I_{D}$.
- `early_voltage_plot()`: this plots $V_{A}$, vs $g_{m}/I_{D}$.

For example, the plot of $I_{D}/W$ vs $g_{m}/I_{D}$ is shown below.

```python
nmos.current_density_plot()
```

![current density plot](./figures/nmos_current_density.svg)

When the lookup table includes a lot of lengths, the plot can become crowded.
You can pass a list of lengths to plot with the `length` parameter.

Use `nmos.lengths` to get a list of all the lengths in the lookup table.

```
array([5.0e-08, 1.0e-07, 2.0e-07, 4.0e-07, 8.0e-07, 1.6e-06, 3.2e-06,
       6.4e-06])
```

Pass a filtered list to the `current_density_plot` method.

```python
nmos.current_density_plot(
    lengths = [5.0e-08, 1.0e-07, 2.0e-07]
)
```

![current density plot](./figures/nmos_current_density_filtered.svg)

Note that the tool does its best to determine how to scale the axes. For
example, in the last plot, a `log` scale was chosen for the y-axis. We can
easily overwrite that, as well as other things.

```python
nmos.current_density_plot(
    lengths = [5.0e-08, 1.0e-07, 2.0e-07],
    y_scale = 'linear',
    x_limit = (5, 20),
    y_limit = (0, 300),
    save_fig="path/to/save/figure/with/extension"
)
```

![current density plot](./figures/nmos_current_density_options.svg)

### Plotting by Expression

Now, suppose we want to plot something completely custom. The example below
shows how.

```python
nmos.plot_by_expression(
    x_axis = nmos.vgs_expression,
    y_axis = {
        "variables": ["id", "gds"],
        "function": lambda x, y: x / y,
        "label": "$I_D / g_{ds} (A/S)$"
        },
)
```

![custom expression](./figures/nmos_custom_expression_1.svg)

For this example, we want $V_{\mathrm{GS}}$ on the x-axis. Since $V_{\mathrm{GS}}$ is such a
commonly-used expression, it is already defined in the code. Other
commonly-used expressions are also defined, such as:

- `gmid_expression`
- `vgs_expression`
- `vds_expression`
- `vsb_expression`
- `gain_expression`
- `current_density_expression`
- `transist_frequency_expression`
- `early_voltage_expression`

For the y-axis, we want a custom expression that uses the parameters $I_D$ and
$g_{\mathrm{ds}}$. This can be done by defining a dictionary that specifies the
variables needed and how to calculate the required parameter. The `label` field
is optional. The function field is also optional if we want to just plot the
parameter, as shown in the example below.

```python
nmos.plot_by_expression(
    x_axis = nmos.vgs_expression,
    # y_axis = nmos.id_expression, ## same as below
    y_axis = {
        "variables": ["id"],
        "label": "$I_D (A)$"
        }
)
```

![custom expression](./figures/nmos_custom_expression_2.svg)

## Looking Up Values

While having plots is a good way to visualize trends, we might also just be
interested in the raw value.

![gain expression](./figures/nmos_gain_plot.svg)

Looking at the figure above, it's hard to read the exact value on the y-axis
for a particular value on the x-axis, especially more so when the scale is
logarithmic. Also, what if we need to read the value for a length that
is not defined in our lookup table?

There are two ways to go about this:

- Zoom in and click on the plot. This prints out the `x` and `y`
  coordinates. Note, in jupyter notebooks, you need to execute `%matplotlib
  widget` or `%matplotlib qt` to interact with the plot.

- Use a lookup method to get a more precise value.

### Lookup Using Interpolation

The snippet below shows how we can lookup some value given the `length` and
`gmid`. The returned value is calculated using interpolation from the
available data. The accuracy of the result depends on how far the points
are from those defined in the table.

```python
x = nmos.lookup_by_gmid(
    length=180e-9,
    gmid=15,
    expression=nmos.gain_expression
)
```

The above code returns a single point. If you want to make a plot of `gain`
vs. `length`, just loop through in the usual manner.

```python
length_sweep = np.arange(100e-9, 1000e-9, 50e-9)
gain_sweep = np.zeros_like(length_sweep)

for idx, length in enumerate(length_sweep):
    gain_sweep[idx] = nmos.lookup_by_gmid(
        length=length,
        gmid=15,
        expression=nmos.gain_expression
    )

nmos.quick_plot(length_sweep, gain_sweep)
```

`lookup_by_gmid()` is a wrapper method over the more flexible method
`lookup_by()`. The example below finds the interpolated length given a
`vgs` and a `gmid`.

```python
x = nmos.lookup_by(
    independent_expression=nmos.vgs_expression,
    independent_value=0.65,
    look_by_expression=nmos.gmid_expression,
    look_by_value=15,
    look_for_expression=nmos.lengths_expression,
)
```

### Lookup By Expression

`lookup_expression_from_table()` simply looks up an expression from the
table. It doesn't use any interpolation. So, make sure that the values you
are looking up are present in the table.

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

The `plot_by_sweep` method is extremely flexible and can be used to create
all sorts of plots. For example, the snippet below shows how to plot the
traditional output characteristic plot of a MOSFET.

```python
nmos.plot_by_sweep(
    lengths=180e-9,
    vsb = 0,
    vds = (0.0, 1, 0.01), # you can also set to `None`
    vgs = (0.0, 1.01, 0.2),
    x_axis_expression = nmos.vds_expression,
    y_axis_expression = nmos.id_expression,
    primary = "vds",
    x_eng_format=True,
    y_eng_format=True,
    y_scale='linear',
)
```

![output characteristic](./figures/nmos_output_characteristics.svg)

### Quick Plot

Let's say we want to see how $V_{\mathrm{DS}_{\mathrm{SAT}}}$ (the drain-source
voltage required to enter saturation) compares with $V_{\mathrm{OV}}$ and
$V^{\star} = \frac{2}{g_m / I_D}$ in a single plot. We can generate each of
these plots individually, as we did before, but ask the method to return the
plot data so that we can combine them in a single plot. Note that you can also
use `lookup_expression_from_table()` to return the required data if you don't
want to see the plot.

```python
vdsat = nmos.plot_by_expression(
    lengths=[45e-9],
    x_axis = nmos.vgs_expression,
    y_axis = nmos.vdsat_expression,
    return_result = True,
)

vov = nmos.plot_by_expression(
    lengths=[45e-9],
    x_axis = nmos.vgs_expression,
    y_axis = {
        "variables": ["vgs", "vth"],
        "function": lambda x, y: x - y,
        },
    return_result = True,
)

vstar = nmos.plot_by_expression(
    lengths=[45e-9],
    x_axis = nmos.vgs_expression,
    y_axis = {
        "variables": ["gm", "id"],
        "function": lambda x, y: 2 / (x/y),
        },
    return_result=True,
)
```

The result is returned in a tuple in the form `(x_data, y_data)`. We can then
make any custom plot using `matplotlib`. Nevertheless, there's a method called
`quick_plot()` that formats the plot in the same way as the generated plots.
`quick_plot()` accepts `numpy` arrays, or a list of `x` and `y` values, as
shown in the example below.

```python
nmos.quick_plot(
    x = [vdsat[0], vstar[0], vov[0]],
    y = [vdsat[1], vstar[1], vov[1]],
    legend = ["$V_{\\mathrm{DS}_{\\mathrm{SAT}}}$", "$V^{\\star}$", "$V_{\\mathrm{OV}}$"],
    x_limit = (0.1, 1),
    y_limit = (0, 0.6),
    x_label = "$V_{\\mathrm{GS}}$",
    y_label = "$V$",
)
```

![qucik plot](./figures/nmos_quick_plot.svg)

# Acknowledgment

- Parsing the output from `hspice` is done using
  [this](https://github.com/HMC-ACE/hspiceParser) script.

- If you find this tool useful, it would be nice if you cite it.
