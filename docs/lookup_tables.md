# Generating a Lookup Table

Everything in Mosplot runs on a lookup table: a pre-computed record of how a
transistor behaves at every bias point. You build one per technology and corner,
save it to a `.npz`, and reuse it for plotting and optimization.

For each model, the generator sweeps a single device and stores every requested
operating-point parameter (`id`, `gm`, `gds`, `cgg`, ...) as a 4-D array indexed
by

```
length × vbs × vgs × vds
```

The device width is fixed (e.g. 10 µm); quantities that scale with width are
normalized by it later.

## 1. Configure a simulator

ngspice, HSPICE, and Spectre take the same arguments; only the class name
changes.

```python
from mosplot.lookup_table_generator.simulators import NgspiceSimulator

sim = NgspiceSimulator(
    simulator_path="ngspice",          # binary name or full path
    temperature=27,                    # °C
    include_paths=[                    # model files to .include ...
        "./NMOS_VTH.inc",
        "./PMOS_VTH.inc",
    ],
    device_parameters={"w": 10e-6},    # fixed width (and any other instance params)
)
```

```python
from mosplot.lookup_table_generator.simulators import HspiceSimulator

sim = HspiceSimulator(
    simulator_path="hspice",
    lib_mappings=[("/path/to/models.lib", "tt")],
    # Transistor wrapped in a subcircuit: (instance name, hierarchical device name)
    mos_spice_symbols=("x1", "x1.main"),
    device_parameters={"w": 10e-6},
)
```

```python
from mosplot.lookup_table_generator.simulators import SpectreSimulator

sim = SpectreSimulator(
    simulator_path="spectre",
    lib_mappings=[("/path/to/models.scs", "tt")],
    device_parameters={"w": 10e-6},
)
```

Arguments:

| Argument             | Default                                | Meaning                                                                              |
| ---                  | ---                                    | ---                                                                                  |
| `simulator_path`     | `"ngspice"` / `"hspice"` / `"spectre"` | Simulator binary. Must be on `PATH` or a full path.                                  |
| `temperature`        | `27`                                   | Simulation temperature in °C.                                                        |
| `include_paths`      | `None`                                 | Model files to include.                                                              |
| `lib_mappings`       | `None`                                 | List of `(library file, section)` pairs, e.g. for process corners.                   |
| `device_parameters`  | `{"w": 10e-6}`                         | Instance parameters written on the transistor line (`w`, `nf`, ...).                 |
| `parameters_to_save` | see below                              | Operating-point parameters to save.                                                  |
| `mos_spice_symbols`  | `("m1", "m1")`                         | Instance name, and the name used to probe it. Change when the model is a subcircuit. |
| `osdi_paths`         | `None`                                 | ngspice only. Compiled Verilog-A models (`.osdi`) to load.                           |
| `hdl_paths`          | `None`                                 | HSPICE only. Verilog-A files to load.                                                |
| `raw_spice`          | `None`                                 | Extra netlist lines, pasted verbatim.                                                |

At least one of `include_paths` or `lib_mappings` is required. Paths are checked
before anything runs.

Available parameters:

| Simulator | Available |
|---|---|
| ngspice | `id`, `weff`, `vth`, `vdsat`, `vdssat`, `gm`, `gmbs`, `gds`, `cgg`, `cgs`, `cbg`, `cgd`, `cdd` |
| HSPICE | `id`, `vth`, `vdsat`, `gm`, `gmbs`, `gds`, `cgg`, `cgs`, `cgd`, `cgb`, `cdd`, `css` |
| Spectre | `id`, `vth`, `vdsat`, `gm`, `gmbs`, `gds`, `cgg`, `cgs`, `cgd`, `cgb`, `cdd`, `css` |

By default everything except `css` is saved. If your model does not expose a
parameter, pass the ones you want in `parameters_to_save` rather than editing
the repo.

## 2. Define the sweeps

One `TransistorSweep` per model. Voltages are `(start, stop, step)` and follow
the device's own sign convention, so PMOS sweeps go negative:

```python
from mosplot.lookup_table_generator import TransistorSweep

nmos_sweep = TransistorSweep(
    mos_type="nmos",
    vgs=(0, 1.2, 0.01),
    vds=(0.01, 1.2, 0.01),
    vbs=(0, -1.2, -0.1),
    length=[130e-9, 180e-9, 250e-9, 500e-9, 1e-6, 2e-6, 5e-6],
)

pmos_sweep = TransistorSweep(
    mos_type="pmos",
    vgs=(0, -1.2, -0.01),
    vds=(-0.01, -1.2, -0.01),
    vbs=(0, 1.2, 0.1),
    length=[130e-9, 180e-9, 250e-9, 500e-9, 1e-6, 2e-6, 5e-6],
)
```

## 3. Build the table

```python
from mosplot.lookup_table_generator import LookupTableGenerator

gen = LookupTableGenerator(
    description="my process, tt, 27C",
    simulator=sim,
    model_sweeps={
        "NMOS_VTH": nmos_sweep,   # key = model name, as it appears in the model file
        "PMOS_VTH": pmos_sweep,
    },
    n_process=1,
)

gen.op_simulation()           # optional: run one operating point and print the netlist
gen.build("./my_process_tt")  # writes ./my_process_tt.npz
```

The keys of `model_sweeps` are the model names written into the netlist and the
names used later to pick a device (`Mosfet(..., mos="NMOS_VTH")`).
`op_simulation()` is a quick sanity check when setting up a new PDK.
`n_process` runs several simulations in parallel; keep it at `1` for ngspice,
which already uses multiple threads internally.

For process corners, build one table per corner, changing only `lib_mappings`
and the file name:

```python
sweeps = {"NMOS_VTH": nmos_sweep, "PMOS_VTH": pmos_sweep}

for corner in ["tt", "ss", "ff"]:
    sim = SpectreSimulator(lib_mappings=[("/path/to/models.scs", corner)])
    gen = LookupTableGenerator(description=corner, simulator=sim, model_sweeps=sweeps)
    gen.build(f"./luts/my_process_{corner}")
```

The optimizer uses these tables to size across all corners at once.

## Loading a table

```python
from mosplot.plot import load_lookup_table

lookup_table = load_lookup_table("./my_process_tt.npz")
print(lookup_table.keys())
# dict_keys(['NMOS_VTH', 'PMOS_VTH', 'description', 'simulator', 'parameter_names', 'device_parameters'])
```

The file is a plain Python dict. Next step: [design charts](plotting.md). A
complete runnable script is in
[`examples/lookup_table_generator`](../examples/lookup_table_generator).
