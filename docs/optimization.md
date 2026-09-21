# Circuit Optimization

The optimizer sizes a complete circuit from its hand-derived operating-point
equations. The circuit is described once, as a Python class. The optimizer then
searches over the design variables, solves the DC operating point of every
candidate using lookups into the [lookup table](lookup_tables.md), evaluates the
performance at every process corner, and writes a netlist of the best design for
verification by simulation.

A library of ready-made designs is included. Each one has a complete circuit
description, so a design in a new technology requires only a short configuration
file with the lookup tables, operating conditions, and target specifications.

- [Using an existing design](#using-an-existing-design)
- [Available designs](#available-designs)
- [Contributing designs](#contributing-designs)
- [Writing a new design](#writing-a-new-design)

## Using an existing design

Every design directory contains two files:


| File               | Contents                                                                | Edited by the user? |
| ---                | ---                                                                     | ---                 |
| `design.py`        | topology, equations, specifications, netlist hooks                      | no                  |
| configuration file | lookup tables, operating conditions, knob bounds, target specifications | yes                 |


`design.py` contains no technology-dependent numbers and is used unchanged. To
target a new technology, copy the example configuration file in the design
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

The adaptation consists of:

1. **Corners:** the lookup table of each corner and the names of the NMOS and
   PMOS devices in it.
2. **Operating conditions:** supply, common-mode and output voltages, load
   capacitance.
3. **Knob bounds:** in particular the minimum and maximum channel length of the
   technology.
4. **Targets:** any specification computed by the design, with a mode
   (`"max"`, `"min"`, or `"eq"`) and a weight.

The run prints the sized transistors and a table of every specification at
every corner, and writes the netlist of the design to `output_module`. The
fields of the configuration file are described in detail in
[Writing a Design from Scratch](writing_designs.md#step-10-the-configuration-file).

## Available designs

In each design, the `nmos` and `pmos` variants refer to the type of the input
pair (amplifiers) or of the mirror and reference devices (current mirrors and
reference circuits).

### Single-ended amplifiers

 | Design                             | Variants                                                                                                                                              |
 | ---                                | ---                                                                                                                                                   |
 | Common source, resistive load      | [NMOS](../examples/designs/VoltageAmplifiers/SingleEnded/CS/nmos/resistive_load)                                                                      |
 | Common source, current-source load | [NMOS](../examples/designs/VoltageAmplifiers/SingleEnded/CS/nmos/pmos_load)                                                                           |
 | Five-transistor OTA                | [NMOS](../examples/designs/VoltageAmplifiers/SingleEnded/5T/nmos) · [PMOS](../examples/designs/VoltageAmplifiers/SingleEnded/5T/pmos)                 |
 | Telescopic cascode OTA             | [NMOS](../examples/designs/VoltageAmplifiers/SingleEnded/Telescopic/nmos) · [PMOS](../examples/designs/VoltageAmplifiers/SingleEnded/Telescopic/pmos) |
 | Folded cascode OTA                 | [NMOS](../examples/designs/VoltageAmplifiers/SingleEnded/Folded/nmos) · [PMOS](../examples/designs/VoltageAmplifiers/SingleEnded/Folded/pmos)         |
 | Rail-to-rail folded cascode OTA    | [complementary input](../examples/designs/VoltageAmplifiers/SingleEnded/RailToRailFolded)                                                             |
 | Two-stage Miller OTA               | [NMOS](../examples/designs/VoltageAmplifiers/SingleEnded/Miller/nmos) · [PMOS](../examples/designs/VoltageAmplifiers/SingleEnded/Miller/pmos)         |

### Fully differential amplifiers


| Design                                        | Variants                                                           |
| ---                                           | ---                                                                |
| Five-transistor OTA with common-mode feedback | [PMOS](../examples/designs/VoltageAmplifiers/Differential/5T/pmos) |


### Current mirrors


 | Design                    | Variants                                                                                                                                    |
 | ---                       | ---                                                                                                                                         |
 | Simple mirror             | [NMOS](../examples/designs/CurrentMirrors/SimpleMirror/nmos) · [PMOS](../examples/designs/CurrentMirrors/SimpleMirror/pmos)                 |
 | Cascode mirror            | [NMOS](../examples/designs/CurrentMirrors/CascodeMirror/nmos) · [PMOS](../examples/designs/CurrentMirrors/CascodeMirror/pmos)               |
 | Wide-swing cascode mirror | [NMOS](../examples/designs/CurrentMirrors/WideSwingMirror/nmos) · [PMOS](../examples/designs/CurrentMirrors/WideSwingMirror/pmos)           |
 | Improved Wilson mirror    | [NMOS](../examples/designs/CurrentMirrors/ImprovedWilsonMirror/nmos) · [PMOS](../examples/designs/CurrentMirrors/ImprovedWilsonMirror/pmos) |


### Reference circuits


| Design                         | Variants                                                                                                                                            |
| ---                            | ---                                                                                                                                                 |
| Constant-$g_m$ bias            | [NMOS](../examples/designs/ReferenceCircuits/ConstantGm/simple/nmos) · [PMOS](../examples/designs/ReferenceCircuits/ConstantGm/simple/pmos)         |
| Cascode constant-$g_m$ bias    | [NMOS](../examples/designs/ReferenceCircuits/ConstantGm/cascode/nmos) · [PMOS](../examples/designs/ReferenceCircuits/ConstantGm/cascode/pmos)       |
| Wide-swing constant-$g_m$ bias | [NMOS](../examples/designs/ReferenceCircuits/ConstantGm/wide-swing/nmos) · [PMOS](../examples/designs/ReferenceCircuits/ConstantGm/wide-swing/pmos) |


The exact specification names of a design are the keys returned by `specs()` in
its `design.py`; all of them appear in the report.

## Contributing designs

Designs for other commonly used analog building blocks are welcome. New designs
can be contributed through a pull request.

## Writing a new design

A topology that is not in the library is described in a new `design.py`: its
connectivity, the variables chosen by the optimizer, the operating-point
equations, and the specifications. The procedure is explained step by step,
using the two-stage Miller OTA as an example, in **[Writing a Design from
Scratch](writing_designs.md)**.
