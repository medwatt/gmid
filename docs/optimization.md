# Circuit Optimization

The optimizer sizes a complete circuit from its hand-derived operating-point
equations. The circuit is described once as a Python class. The optimizer then
searches over the design variables, solves the DC operating point of every
candidate using lookups into the [lookup table](lookup_tables.md), evaluates the
performance at every process corner, and writes a netlist of the best design for
verification by simulation.

A library of ready-made designs is included, so a design in a new technology
needs only a short configuration file with the tables, operating conditions, and
targets.

## Using an existing design

Every design directory contains two files: `design.py`, which holds the
topology, equations, specifications, and is never edited; and a configuration
file with the lookup tables, operating conditions, knob bounds, and targets,
which the user adapts.

To target a new technology, copy the example configuration file in the design
directory, adapt it, and run it:

```python
from design import Circuit, run
from mosplot.optimizer import Corner, Knob, Spec

# operating conditions read by the circuit model
COND = dict(vdd=1.2, vin_cm=0.4, vout_dc=0.6, cout=5e-12)
LMIN, LMAX = 130e-9, 2e-6

# search range of every knob declared in design.py
PARAMETERS = [
    Knob("M1a_GMID", (10, 20)),
    Knob("M1a_L", (LMIN, LMAX)),
    Knob("M1a_ID", (5e-6, 30e-6)),
    # ...
]

# targets; the keys are the specifications computed by design.py
TARGET_SPECS = {
    "GBW":          Spec(10e6,  "max", 1.0),
    "AC Gain (dB)": Spec(50.0,  "max", 1.0),
    "PM":           Spec(70.0,  "max", 2.0),
    "Itotal":       Spec(50e-6, "min", 1.0),
}

# one lookup table per process corner; the first corner is the reference
CORNERS = [
    Corner("tt", "luts/tt.npz", "nmos", "pmos", conditions=COND),
    Corner("ss", "luts/ss.npz", "nmos", "pmos", conditions=COND),
    Corner("ff", "luts/ff.npz", "nmos", "pmos", conditions=COND),
]

if __name__ == "__main__":
    run(
        circuit=Circuit,
        parameters=PARAMETERS,
        target_specs=TARGET_SPECS,
        corners=CORNERS,
        output_module="./design.scs",
        workers=4,
    )
```

The adaptation consists of the corners (table and device names per corner), the
operating conditions (supply, common-mode and output voltages, load capacitance),
the knob bounds (in particular the technology's minimum and maximum length), and
the targets (any spec computed by the design, with a mode and a weight).

The run prints the sized transistors and a table of every spec at every corner,
and writes the netlist to `output_module`. The configuration fields are described
in [writing a design](writing_designs.md).

## Available designs

In each design the `nmos` and `pmos` variants refer to the type of the input
pair (amplifiers) or of the mirror and reference devices (current mirrors and
reference circuits).

Single-ended amplifiers:

| Design | Variants |
|---|---|
| Common source, resistive load | [NMOS](../examples/designs/VoltageAmplifiers/SingleEnded/CS/nmos/resistive_load) |
| Common source, current-source load | [NMOS](../examples/designs/VoltageAmplifiers/SingleEnded/CS/nmos/pmos_load) |
| Five-transistor OTA | [NMOS](../examples/designs/VoltageAmplifiers/SingleEnded/5T/nmos) · [PMOS](../examples/designs/VoltageAmplifiers/SingleEnded/5T/pmos) |
| Telescopic cascode OTA | [NMOS](../examples/designs/VoltageAmplifiers/SingleEnded/Telescopic/nmos) · [PMOS](../examples/designs/VoltageAmplifiers/SingleEnded/Telescopic/pmos) |
| Folded cascode OTA | [NMOS](../examples/designs/VoltageAmplifiers/SingleEnded/Folded/nmos) · [PMOS](../examples/designs/VoltageAmplifiers/SingleEnded/Folded/pmos) |
| Rail-to-rail folded cascode OTA | [complementary input](../examples/designs/VoltageAmplifiers/SingleEnded/RailToRailFolded) |
| Two-stage Miller OTA | [NMOS](../examples/designs/VoltageAmplifiers/SingleEnded/Miller/nmos) · [PMOS](../examples/designs/VoltageAmplifiers/SingleEnded/Miller/pmos) |

Fully differential amplifiers:

| Design | Variants |
|---|---|
| Five-transistor OTA with common-mode feedback | [PMOS](../examples/designs/VoltageAmplifiers/Differential/5T/pmos) |

Current mirrors:

| Design | Variants |
|---|---|
| Simple mirror | [NMOS](../examples/designs/CurrentMirrors/SimpleMirror/nmos) · [PMOS](../examples/designs/CurrentMirrors/SimpleMirror/pmos) |
| Cascode mirror | [NMOS](../examples/designs/CurrentMirrors/CascodeMirror/nmos) · [PMOS](../examples/designs/CurrentMirrors/CascodeMirror/pmos) |
| Wide-swing cascode mirror | [NMOS](../examples/designs/CurrentMirrors/WideSwingMirror/nmos) · [PMOS](../examples/designs/CurrentMirrors/WideSwingMirror/pmos) |
| Improved Wilson mirror | [NMOS](../examples/designs/CurrentMirrors/ImprovedWilsonMirror/nmos) · [PMOS](../examples/designs/CurrentMirrors/ImprovedWilsonMirror/pmos) |

Reference circuits:

| Design | Variants |
|---|---|
| Constant-gm bias | [NMOS](../examples/designs/ReferenceCircuits/ConstantGm/simple/nmos) · [PMOS](../examples/designs/ReferenceCircuits/ConstantGm/simple/pmos) |
| Cascode constant-gm bias | [NMOS](../examples/designs/ReferenceCircuits/ConstantGm/cascode/nmos) · [PMOS](../examples/designs/ReferenceCircuits/ConstantGm/cascode/pmos) |
| Wide-swing constant-gm bias | [NMOS](../examples/designs/ReferenceCircuits/ConstantGm/wide-swing/nmos) · [PMOS](../examples/designs/ReferenceCircuits/ConstantGm/wide-swing/pmos) |

The exact spec names of a design are the keys returned by `specs()` in its
`design.py`; all of them appear in the report.

## Contributing designs

Designs for other commonly used analog building blocks are welcome through a
pull request.

## Writing a new design

A topology that is not in the library is described in a new `design.py`: its
connectivity, the variables chosen by the optimizer, the operating-point
equations, and the specifications. The procedure is explained step by step,
using the two-stage Miller OTA as an example, in [Writing a Design from
Scratch](writing_designs.md).
