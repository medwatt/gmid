# Writing a Design from Scratch

This page explains how the optimizer works and builds a complete design file
from scratch: a two-stage Miller OTA with a PMOS input pair. The same procedure
applies to any other topology. The finished files are in
[`examples/designs/VoltageAmplifiers/SingleEnded/Miller/pmos`](../examples/designs/VoltageAmplifiers/SingleEnded/Miller/pmos).

To optimize one of the [existing designs](optimization.md#available-designs),
none of this is needed; only a configuration file is written.

- [How the optimizer works](#how-the-optimizer-works)
- [Worked example: the Miller OTA](#worked-example-the-miller-ota)
  - [Step 1: Topology](#step-1-topology)
  - [Step 2: Knobs](#step-2-knobs)
  - [Step 3: Derive the node voltages](#step-3-derive-the-node-voltages)
  - [Step 4: Unknowns](#step-4-unknowns)
  - [Step 5: `solve_point()`](#step-5-solve_point)
  - [Step 6: `residuals()`](#step-6-residuals)
  - [Step 7: `specs()`](#step-7-specs)
  - [Step 8: Multicorner conservation](#step-8-multicorner-conservation)
  - [Step 9: Netlist hooks](#step-9-netlist-hooks)
  - [Step 10: The configuration file](#step-10-the-configuration-file)
  - [Step 11: Run and read the results](#step-11-run-and-read-the-results)
- [Checklist for a new circuit](#checklist-for-a-new-circuit)
- [Parallel evaluation](#parallel-evaluation)

## How the optimizer works

### Three kinds of variables

Every quantity in a circuit model belongs to one of three groups.

 | Group                  | Chosen by                    | Examples                                                                                              |
 | ---                    | ---                          | ---                                                                                                   |
 | **Knobs**              | the outer optimizer (CMA-ES) | $\mathrm{gm}/I_D$, channel length, a branch current, a mirror ratio, $C_C$, $R_Z$                     |
 | **Unknowns**           | an inner nonlinear solver    | a node voltage that depends on itself through a lookup, e.g. the $V_{DS}$ of a diode-connected device |
 | **Derived quantities** | computed directly            | widths, currents, node voltages that can be written explicitly                                        |

Widths are never knobs. For a device biased at a given $\mathrm{gm}/I_D$,
length, and $V_{DS}$, the lookup table returns the current density $I_D/W$, and
the width follows as $W = I_D / (I_D/W)$. This is the $\mathrm{gm}/I_D$ method
applied automatically to every transistor.

### One evaluation of a candidate design

```
knobs ──> size at the reference corner ──> freeze the geometry ──> re-solve every other corner ──> worst case of each spec ──> cost
            (solve the unknowns,              (W, L, passives,        (same W and L; gm/Id and
             derive W, compute specs)          conserved biases)       currents move with process)
```

1. **Sizing.** At the reference corner (the first in the list), the knobs are
   fixed, the unknowns are solved, and the widths are derived.
2. **Freezing.** Widths, lengths, passive values, and conserved bias quantities
   are stored. These define the physical circuit.
3. **Re-cornering.** At every other corner, the physical circuit stays fixed and
   its operating point is solved again. A transistor of frozen width and length
   runs at a different $\mathrm{gm}/I_D$ in the fast corner than in the slow one; the solver
   finds that $\mathrm{gm}/I_D$.
4. **Cost.** For each spec, the worst value over all corners is compared with
   its target.

The outer search uses CMA-ES followed by an SLSQP polish. Candidates whose DC
operating point cannot be solved receive a large penalty, graded by how far the
solver was from convergence.

### Two files per design

| File               | Contents                                                        | Depends on the technology? |
| ---                | ---                                                             | ---                        |
| `design.py`        | topology, knob names, unknowns, equations, specs, netlist hooks | no                         |
| configuration file | lookup tables, operating conditions, knob bounds, target specs  | yes                        |

The same `design.py` is reused unchanged in any technology; only the
configuration file changes.

## Worked example: the Miller OTA

![Miller OTA](figures/miller_ota.svg)

The circuit is a two-stage OTA:

- **First stage:** PMOS input pair `M1a`/`M1b`, NMOS current-mirror load
  `M2a`/`M2b`, PMOS tail current source `M3`.
- **Second stage:** NMOS common-source device `M4` with PMOS current-source load
  `M5`.
- **Compensation:** Miller capacitor `CC` with nulling resistor `Rz`.
- **Bias:** the gates of `M3` and `M5` share the bias node `vbp`. In the
  generated netlist, `vbp` is produced by a diode-connected replica `Mvbp` fed
  by a reference current, exactly as on silicon.

The colors in the schematic map directly onto the code: blue names are knobs,
green names are unknowns, and red expressions are node voltages computed in
`solve_point()`.

The file starts with the imports:

```python
from __future__ import annotations

import numpy as np

from mosplot.optimizer import (
    CircuitModel, Instance, Knob, Passive, Spec, State, Unknown, VSource,
    build_ss_model, run, vres, rres,
)
```

### Step 1: Topology

The topology is a direct transcription of the schematic shown above.

```python
class Circuit(CircuitModel):
    NAME = "amp"                                    # subcircuit name in the netlist
    PORTS = ["VINN", "VINP", "VOUT", "vdd", "vss"]  # subcircuit port order
    GROUND = "vss"

    MOSFETS = [
        Instance("M1a", "pmos", d="n01",  g="VINN", s="n02", b="vdd"),
        Instance("M1b", "pmos", d="n03",  g="VINP", s="n02", b="vdd"),
        Instance("M2a", "nmos", d="n01",  g="n01",  s="vss", b="vss"),
        Instance("M2b", "nmos", d="n03",  g="n01",  s="vss", b="vss"),
        Instance("M3",  "pmos", d="n02",  g="vbp",  s="vdd", b="vdd"),
        Instance("M4",  "nmos", d="VOUT", g="n03",  s="vss", b="vss"),
        Instance("M5",  "pmos", d="VOUT", g="vbp",  s="vdd", b="vdd"),
    ]
    PASSIVES = [
        Passive("Rz",   "res", a="n04",  b="n03"),
        Passive("CC",   "cap", a="n04",  b="VOUT"),
        Passive("COUT", "cap", a="VOUT", b="vss", external=True),
    ]
    VSOURCES = [
        VSource("VDD", p="vdd", n="vss", supply=True),
        VSource("VBP", p="vbp", n="vss", mirror="M3"),
    ]
    SIGNAL_NODES = {"VINP", "VINN"}
```


| Declaration                                 | Meaning                                                                                                                                                                                                                                                |
| ---                                         | ---                                                                                                                                                                                                                                                    |
| `Instance(name, "nmos"/"pmos", d, g, s, b)` | One transistor and its four terminals.                                                                                                                                                                                                                 |
| `Passive(name, "res"/"cap", a, b)`          | A resistor or capacitor. `external=True` marks an element that belongs to the small-signal model but not to the circuit, such as the load capacitance; its value comes from the operating conditions.                                                  |
| `VSource(name, p, n, supply=True)`          | The supply.                                                                                                                                                                                                                                            |
| `VSource(name, p, n, mirror="M3")`          | A bias node realized by a current mirror. In the small-signal model the node is an AC ground. In the netlist, it is written as a diode-connected copy of `M3` fed by a reference current. The replica `Mvbp` is therefore **not** listed in `MOSFETS`. |
| `SIGNAL_NODES`                              | Nodes driven by the AC input. All other nodes attached to a voltage source are AC grounds.                                                                                                                                                             |

Node names are free-form, but every net tied to the lower rail must use the
name given in `GROUND`.

### Step 2: Knobs

Knobs are declared by name and **role** only. Their numeric bounds belong to
the configuration file, which keeps `design.py` free of technology numbers.

```python
    KNOBS = [
        Knob("M1a_GMID", role="op", sets_width_of="M1a"),
        Knob("M2a_GMID", role="op", sets_width_of="M2a"),
        Knob("M3_GMID",  role="op", sets_width_of="M3"),
        Knob("M4_GMID",  role="op", sets_width_of="M4"),
        Knob("M1a_L", role="geom"),
        Knob("M2a_L", role="geom"),
        Knob("M3_L",  role="geom"),
        Knob("M4_L",  role="geom"),
        Knob("W5_over_W3", role="geom"),
        Knob("CC", role="geom"),
        Knob("Rz", role="geom"),
        Knob("M1a_ID", role="external"),
    ]
```

The role determines what happens to the knob at corners other than the
reference:

| Role         | Use for                                                      | At the reference corner | At other corners                                                                             |
| ---          | ---                                                          | ---                     | ---                                                                                          |
| `"op"`       | $\mathrm{gm}/I_D$ of a device                                | chosen by the optimizer | re-solved so that the device width `sets_width_of` equals its frozen value                   |
| `"geom"`     | lengths, width ratios, multipliers, passive values           | chosen by the optimizer | held fixed                                                                                   |
| `"external"` | externally applied bias: a reference current, a bias voltage | chosen by the optimizer | held fixed, unless listed in `RECORNER_RESOLVE` ([Step 8](#step-8-multicorner-conservation)) |

Some guidelines:

- **Matched devices share knobs.** `M1a` and `M1b` form a differential pair, so
  both use `M1a_GMID` and `M1a_L`. Likewise `M2a`/`M2b`.
- **Mirror copies share $\mathrm{gm}/I_D$ and length with their master.** `M5` copies `M3`,
  so it uses `M3_GMID` and `M3_L`; its size is set by the ratio `W5_over_W3`.
- **One current knob per independent branch.** Here `M1a_ID` sets the
  first-stage branch current; every other current follows from KCL and mirror
  ratios.

The knob names are arbitrary. They are read back as attributes, e.g.
`v.M1a_GMID`, and must match the names in the configuration file.

### Step 3: Derive the node voltages

This step is done on paper. Every transistor needs a $V_{DS}$ (and a $V_{SB}$)
for its lookup. The lookups take and return magnitudes, so every voltage in the
model is positive. For a PMOS device, $V_{DS}$ and $V_{GS}$ are negative, and
the values used are $|V_{DS}| = V_S - V_D$ and $|V_{GS}| = V_S - V_G$. The
lookups apply the signs internally.

Start from the nodes whose voltages are known:

| Node   | Voltage                  | Reason                                                      |
| ---    | ---                      | ---                                                         |
| `VOUT` | $V_{out,dc}$             | the intended output operating point, an operating condition |
| `n03`  | $V_{GS,M4}$              | gate of `M4`, whose source is at ground                     |
| `n02`  | $V_{in,cm} + V_{GS,M1b}$ | source of the input pair, one $V_{GS}$ above the input      |
| `n01`  | $V_{GS,M2}$              | gate of the mirror, whose source is at ground               |
| `vbp`  | $V_{DD} - V_{GS,M3}$     | gate of the tail source                                     |

Then write each $V_{DS}$ and check whether its right-hand side can be computed
**before** the device itself is looked up:

| Device | $V_{DS}$                              | Computable in advance?                                |
| ---    | ---                                   | ---                                                   |
| `M4`   | $V_{out,dc}$                          | yes: a condition                                      |
| `M5`   | $V_{DD} - V_{out,dc}$                 | yes: a condition                                      |
| `M2b`  | $V_{GS,M4}$                           | yes, once `M4` has been looked up                     |
| `M2a`  | $V_{GS,M2a}$ (diode)                  | **no**: it needs its own $V_{GS}$ → unknown `M2a_VDS` |
| `M1b`  | $V_{in,cm} + V_{GS,M1b} - V_{GS,M4}$  | **no**: it needs its own $V_{GS}$ → unknown `M1b_VDS` |
| `M1a`  | $V_{in,cm} + V_{GS,M1b} - V_{GS,M2b}$ | yes, once `M1b` and `M2b` are known                   |
| `M3`   | $V_{DD} - V_{in,cm} - V_{GS,M1b}$     | yes, once `M1b` is known                              |


A self-referential equation such as $V_{DS,M2a} = V_{GS,M2a}(V_{DS,M2a})$ cannot
be evaluated in a single pass, so it becomes an unknown. Everything else is
evaluated directly, in the right order. This analysis fixes both the unknowns
and the order of the lookups in `solve_point()`.

### Step 4: Unknowns

```python
    UNKNOWNS = [
        Unknown("M1b_VDS",   seed=lambda c: c["vdd"] / 3, bound=lambda c: (0.02, c["vdd"])),
        Unknown("M2a_VDS",   seed=lambda c: c["vdd"] / 3, bound=lambda c: (0.02, c["vdd"])),
        Unknown("Mvbp_GMID", seed=lambda c: 12.0,         bound=lambda c: (4.0, 30.0)),
    ]
```

Each unknown has an initial guess (`seed`) and a search box (`bound`). Both can
be constants or functions of the operating conditions `c`. Writing them in
terms of `c["vdd"]` keeps the design independent of the supply voltage.

The third unknown belongs to the bias replica. `Mvbp` has the same width and
length as `M3`, is diode-connected, and must produce the same gate voltage as
`M3`. Its $V_{DS}$ equals its own $V_{GS}$, whereas `M3` operates at a different
$V_{DS}$, so the two devices do not run at the same $\mathrm{gm}/I_D$. The solver finds the
$\mathrm{gm}/I_D$ at which the replica's $V_{GS}$ equals that of `M3`. From it follows the
reference current that the replica needs, which is written into the netlist.

> [!NOTE]
> Keep the number of unknowns small. A variable that can be computed directly
> should not be an unknown: it slows the solver down and can make it less
> robust.

### Step 5: `solve_point()`

`solve_point()` performs every lookup and derives every quantity. It is called
many times per candidate by the inner solver, each time with new trial values
of the unknowns.

```python
    def solve_point(self, v: State, dev, cond) -> State:
```

| Argument | Contents                                                                                               |
| ---      | ---                                                                                                    |
| `v`      | the knobs and the current trial values of the unknowns, as attributes (`v.M1a_GMID`, `v.M2a_VDS`, ...) |
| `dev`    | the lookup tables of the active corner: `dev.nmos(...)` and `dev.pmos(...)`                            |
| `cond`   | the operating conditions of the active corner (`cond["vdd"]`, ...)                                     |

**Operating conditions.**

```python
COUT = cond["cout"]
VDD = cond["vdd"]
VIN_CM = cond["vin_cm"]
VOUT_DC = cond["vout_dc"]
```

**Lookups, in the order found in Step 3.**

```python
M4_VDS = VOUT_DC
M4 = dev.nmos(gmid=v.M4_GMID, L=v.M4_L, vds=M4_VDS, vsb=0.0)
M2a = dev.nmos(gmid=v.M2a_GMID, L=v.M2a_L, vds=v.M2a_VDS, vsb=0.0)   # unknown VDS
M2b_VDS = M4.vgs
M2b = dev.nmos(gmid=v.M2a_GMID, L=v.M2a_L, vds=M2b_VDS, vsb=0.0)
M1b = dev.pmos(gmid=v.M1a_GMID, L=v.M1a_L, vds=v.M1b_VDS, vsb=0.0)   # unknown VDS
M1a_VDS = VIN_CM + M1b.vgs - M2b.vgs
M1a = dev.pmos(gmid=v.M1a_GMID, L=v.M1a_L, vds=M1a_VDS, vsb=0.0)
M3_VDS = VDD - VIN_CM - M1b.vgs
M3 = dev.pmos(gmid=v.M3_GMID, L=v.M3_L, vds=M3_VDS, vsb=0.0)
M5_VDS = VDD - VOUT_DC
M5 = dev.pmos(gmid=v.M3_GMID, L=v.M3_L, vds=M5_VDS, vsb=0.0)
Mvbp = dev.pmos(gmid=v.Mvbp_GMID, L=v.M3_L, vds=M3.vgs, vsb=0.0)    # diode: VDS = VGS of M3
```

Each call returns a `DevicePoint`. All of its fields are positive magnitudes:

| Field               | Quantity                         |
| ---                 | ---                              |
| `vgs`               | $\lvert V_{GS} \rvert$           |
| `vdsat`             | $\lvert V_{DS,sat} \rvert$       |
| `jd`                | current density $I_D/W$ (A/m)    |
| `gds_id`            | $g_{ds}/I_D$                     |
| `cgs`, `cgd`, `cdd` | capacitances of the table device |
| `vds_used`          | the $V_{DS}$ used for the lookup |


`M2a` and `M2b` share a gate voltage, but they are looked up separately because
their drain voltages differ, and the lookup table captures that difference.

**Currents** follow from KCL and the mirror ratios:

```python
ID = {}
ID["M1a"] = v.M1a_ID
ID["M1b"] = ID["M1a"]
ID["M2a"] = ID["M1a"]
ID["M2b"] = ID["M2a"]
ID["M3"] = ID["M1a"] + ID["M1b"]
ID["M5"] = ID["M3"] * v.W5_over_W3
ID["M4"] = ID["M5"]
ID["Mvbp"] = ID["M3"]
IREF_Mvbp = ID["M3"] * Mvbp.jd / M3.jd
```

`IREF_Mvbp` is the current of a diode with the width of `M3` at the gate voltage
of `M3`: $W_3 \cdot (I_D/W)_{Mvbp}$. It is the reference current the bias
network must supply.

**Widths** come from the current density, and matched devices copy the width of
their partner:

```python
W = {}
W["M1a"] = ID["M1a"] / M1a.jd
W["M1b"] = W["M1a"]
W["M2a"] = ID["M2a"] / M2a.jd
W["M2b"] = W["M2a"]
W["M3"] = ID["M3"] / M3.jd
W["M4"] = ID["M4"] / M4.jd
W["M5"] = W["M3"] * v.W5_over_W3
W["Mvbp"] = W["M3"]
```

**Lengths and $\mathrm{gm}/I_D$** for the report and the netlist:

```python
L = {"M1a": v.M1a_L, "M1b": v.M1a_L, "M2a": v.M2a_L, "M2b": v.M2a_L,
     "M3": v.M3_L, "M4": v.M4_L, "M5": v.M3_L, "Mvbp": v.M3_L}
GMID = {"M1a": v.M1a_GMID, "M1b": v.M1a_GMID, "M2a": v.M2a_GMID, "M2b": v.M2a_GMID,
        "M3": v.M3_GMID, "M4": v.M4_GMID, "M5": v.M3_GMID, "Mvbp": v.Mvbp_GMID}
```

**Small-signal parameters** for every device in `MOSFETS`. `small_signal()`
scales the table capacitances from the width of the table device to the actual
width:

```python
ptw, ntw = dev.pmos.table_width, dev.nmos.table_width
ss = {
    "M1a": M1a.small_signal(v.M1a_GMID, ID["M1a"], W["M1a"], ptw, use_gmb=False),
    "M1b": M1b.small_signal(v.M1a_GMID, ID["M1b"], W["M1b"], ptw, use_gmb=False),
    "M2a": M2a.small_signal(v.M2a_GMID, ID["M2a"], W["M2a"], ntw, use_gmb=False),
    "M2b": M2b.small_signal(v.M2a_GMID, ID["M2b"], W["M2b"], ntw, use_gmb=False),
    "M3":  M3.small_signal(v.M3_GMID,  ID["M3"],  W["M3"],  ptw, use_gmb=False),
    "M4":  M4.small_signal(v.M4_GMID,  ID["M4"],  W["M4"],  ntw, use_gmb=False),
    "M5":  M5.small_signal(v.M3_GMID,  ID["M5"],  W["M5"],  ptw, use_gmb=False),
}
```

`use_gmb=True` includes the body transconductance, which matters for cascodes.
All sources in this circuit are at their bulk, so it is disabled.

**Return** every quantity that `residuals()`, `specs()`, or the netlist hooks
will read. `W`, `L`, `ID`, `GMID`, and `ss` are required by the framework:

```python
return State(
    W=W, L=L, ID=ID, GMID=GMID, ss=ss,
    M1a=M1a, M1b=M1b, M2a=M2a, M2b=M2b, M3=M3, M4=M4, M5=M5, Mvbp=Mvbp,
    IREF_Mvbp=IREF_Mvbp,
    VDD=VDD, VIN_CM=VIN_CM, VOUT_DC=VOUT_DC, COUT=COUT, VBP=VDD - M3.vgs,
    CC=v.CC, Rz=v.Rz,
    M1b_VDS=v.M1b_VDS, M2a_VDS=v.M2a_VDS, Mvbp_GMID=v.Mvbp_GMID,
)
```

### Step 6: `residuals()`

One equation per unknown, each written as "left side minus right side":

```python
def residuals(self, b) -> list:
    return [
        vres(b.M1b_VDS, b.VIN_CM + b.M1b.vgs - b.M4.vgs),   # V(n02) - V(n03)
        vres(b.M2a_VDS, b.M2a.vgs),                         # diode: VDS = VGS
        vres(b.Mvbp.vgs, b.M3.vgs),                         # replica matches M3
    ]
```

`b` is the `State` returned by `solve_point()`. The helpers normalize the
residuals so that different equations carry comparable weight:

| Helper                      | Value                           | Use for                               |
| ---                         | ---                             | ---                                   |
| `vres(lhs, rhs, scale=0.5)` | $(lhs - rhs)/0.5\,\mathrm{V}$   | voltage equalities                    |
| `rres(lhs, rhs, ref)`       | $(lhs - rhs)/\lvert ref \rvert$ | relative equalities, such as currents |

The number of residuals must equal the number of unknowns. The optimizer
checks this when it starts and raises an error otherwise.

### Step 7: `specs()`

`specs()` turns the solved state into performance figures. It is a pure
function of `b`: no lookups and no solving.

**AC performance.** `build_ss_model()` assembles the small-signal circuit from
the topology, the per-device parameters, and the passive values. `transfer()`
then applies a stimulus to the input nodes and observes a weighted sum of output
nodes:

```python
def specs(self, b, cond) -> dict:
    ss = build_ss_model(
        self.MOSFETS, self.PASSIVES, self.VSOURCES, b.ss,
        {"Rz": b.Rz, "CC": b.CC, "COUT": cond["cout"]},
        signal_nodes=self.SIGNAL_NODES,
    )
    ac = ss.transfer(inputs={"VINP": 0.5, "VINN": -0.5}, output={"VOUT": 1.0})   # differential
    ac_cm = ss.transfer(inputs={"VINP": 1.0, "VINN": 1.0}, output={"VOUT": 1.0})  # common mode
    out = {
        "GBW": ac.ugf(),
        "AC Gain (dB)": 20.0 * np.log10(max(abs(ac.gain()), 1e-300)),
        "PM": ac.phase_margin(),
        "DC CMR (dB)": 20.0 * np.log10(max(ac.rejection(ac_cm), 1e-300)),
    }
```

| Method             | Returns                                              |
| ---                | ---                                                  |
| `gain()`           | DC gain                                              |
| `ugf()`            | unity-gain frequency (Hz)                            |
| `phase_margin()`   | phase margin (degrees)                               |
| `response(f)`      | complex response at frequency `f`                    |
| `rejection(other)` | ratio of this gain to the gain of `other`, e.g. CMRR |

A fully differential output is observed with `output={"VOUTP": 1.0, "VOUTN": -1.0}`.

**Output centering.** The output voltage $V_{out,dc}$ is an input of the
model: it sets the $V_{DS}$ of `M4` and `M5`. Nothing yet guarantees that the
first stage actually delivers the gate voltage `M4` needs for that output. In
the real circuit, `n03` sits at the mirror voltage $V_{GS,M2}$, while `M4`
requires $V_{GS,M4}$. The difference is the systematic offset of the amplifier.
If it is not small, the output of the open-loop circuit saturates at a rail in
simulation:

```python
out["VOUT_Error"] = abs(b.M4.vgs - b.M2a.vgs)
```

This quantity is given a high weight in the configuration file.

**Area, current, bias voltage, and signal ranges.** Every saturation condition
$V_{DS} \geq V_{DS,sat}$ turns into a limit on the input or output voltage:

```python
out["Area"] = sum(b.L[m] * b.W[m] for m in b.W)
out["Itotal"] = b.ID["M1a"] + b.ID["M1b"] + b.ID["M5"]     # current drawn from VDD
out["VBP"] = b.VDD - b.M3.vgs
out["VOUT_MAX"] = b.VDD - b.M5.vdsat
out["VOUT_MIN"] = b.M4.vdsat
out["VIN_MAX"] = b.VDD - b.M3.vdsat - b.M1b.vgs
out["VIN_MIN"] = max(
    b.M1a.vdsat - b.M1b.vgs + b.M2b.vgs,
    b.M1b.vdsat - b.M1b.vgs + b.M4.vgs,
)
out["Output_Swing"] = out["VOUT_MAX"] - out["VOUT_MIN"]
return out
```

Every key returned here appears in the report. Only the keys listed in the
configuration file's `TARGET_SPECS` affect the optimization.

### Step 8: Multicorner conservation

On silicon, the bias of this amplifier is a fixed reference current flowing
into the diode `Mvbp`. When the process changes, that current stays the same,
while `vbp` and the tail current of `M3` move. The model must behave in the
same way at every non-reference corner.

Three declarations express this:

```python
RECORNER_RESOLVE = ["M1a_ID"]

def freeze_extra(self, b) -> dict:
    return {"IREF_Mvbp": b.IREF_Mvbp}

def recorner_residuals(self, b, frozen) -> list:
    e = frozen["extra"]
    return [rres(b.IREF_Mvbp, e["IREF_Mvbp"], e["IREF_Mvbp"])]
```

1. `freeze_extra()` stores the reference current computed at the reference
   corner.
2. At the other corners, `RECORNER_RESOLVE` releases the knob `M1a_ID`, which
   is normally held fixed because of its `"external"` role.
3. `recorner_residuals()` adds the condition that pins it: the reference
   current must equal its frozen value. The branch current therefore drifts
   with the process, exactly as in the real mirror.

The re-solve must remain square: each entry in `recorner_residuals()` needs
exactly one knob in `RECORNER_RESOLVE`. The optimizer checks this when it
starts.

The same pattern applies to a bias applied as a fixed voltage (a `VSource`
without `mirror`) whose value is derived from a knob. Store the voltage in
`freeze_extra()`, add a `vres(...)` for it in `recorner_residuals()`, and
release the knob that determines it. A circuit without derived biases does not
need these hooks.

### Step 9: Netlist hooks

After optimization, the frozen design is written as a Spectre subcircuit. The
hooks provide the values that the topology alone does not contain:

```python
def mirror_currents(self, ref_op) -> dict:
    return {"VBP": ref_op.IREF_Mvbp}          # reference current of the vbp replica

def passive_values(self, ref_op) -> dict:
    return {"Rz": ref_op.Rz, "CC": ref_op.CC}  # optimized passives

def netlist_context(self, corner, ref_op=None) -> dict:
    return {"vcm": corner.cond("vin_cm"), "vout_dc": corner.cond("vout_dc")}
```

| Hook                             | Purpose                                                                                          |
| ---                              | ---                                                                                              |
| `passive_values()`               | values of every non-external `Passive`                                                           |
| `mirror_currents()`              | reference current of every `mirror=` bias. Without it, the current of the master device is used. |
| `vsource_values(ref_op, frozen)` | DC value of every bias written as an ideal voltage source                                        |
| `isource_values(ref_op, frozen)` | DC value of every entry in `ISOURCES`                                                            |
| `netlist_context()`              | values written as netlist `parameters`, for use by a testbench                                   |
| `extra_netlist_lines()`          | lines copied verbatim into the subcircuit                                                        |

### Step 10: The configuration file

The configuration file binds the design to a technology. It supplies the
lookup tables, the operating conditions, the knob bounds, and the targets.

```python
from design import Circuit, run
from mosplot.optimizer import Corner, Knob, Spec

COND = dict(vdd=1.2, vin_cm=0.4, vout_dc=0.6, cout=5e-12)
LMIN, LMAX = 130e-9, 2e-6

PARAMETERS = [
    Knob("M1a_GMID", (10, 20)),
    Knob("M2a_GMID", (10, 20)),
    Knob("M3_GMID", (10, 20)),
    Knob("M4_GMID", (10, 20)),
    Knob("M1a_L", (LMIN, LMAX)),
    Knob("M2a_L", (LMIN, LMAX)),
    Knob("M3_L", (LMIN, LMAX)),
    Knob("M4_L", (LMIN, LMAX)),
    Knob("W5_over_W3", (1, 16)),
    Knob("CC", (1e-12, 2e-12)),
    Knob("Rz", (100, 10e3)),
    Knob("M1a_ID", (5e-6, 30e-6)),
]

TARGET_SPECS = {
    "GBW":          Spec(10e6,   "max", 1.0),
    "AC Gain (dB)": Spec(50.0,   "max", 1.0),
    "DC CMR (dB)":  Spec(50.0,   "max", 0.5),
    "PM":           Spec(70.0,   "max", 2.0),
    "Area":         Spec(50e-12, "min", 0.3),
    "Itotal":       Spec(50e-6,  "min", 1.0),
    "VOUT_Error":   Spec(0.02,   "min", 20.0),
}

CORNERS = [
    Corner("tt", "luts/tt.npz", "nmos", "pmos", conditions=COND),
    Corner("ff", "luts/ff.npz", "nmos", "pmos", conditions=COND),
    Corner("ss", "luts/ss.npz", "nmos", "pmos", conditions=COND),
]

if __name__ == "__main__":
    run(
        circuit=Circuit,
        parameters=PARAMETERS,
        target_specs=TARGET_SPECS,
        corners=CORNERS,
        output_module="./design.scs",
        maxiter=100,
        seed=1,
        workers=4,
    )
```

**Operating conditions.** `COND` holds everything `solve_point()` reads from
`cond`. A voltage or temperature corner reuses a table with different
conditions, e.g. `Corner("tt_lowvdd", "luts/tt.npz", "nmos", "pmos", conditions=dict(COND, vdd=1.08))`.

**Knob bounds.** Every knob of `design.py` needs a bound. A tuple `(lo, hi)` is
a continuous range.

**Specs.** `Spec(target, mode, weight)`:

| Mode    | Meaning                             | Penalty                                                                                                                |
| ---     | ---                                 | ---                                                                                                                    |
| `"max"` | larger is better, at least `target` | grows with $\log(\text{target}/\text{value})$ below the target                                                         |
| `"min"` | smaller is better, at most `target` | grows with the relative excess above the target; a small term keeps pushing the value down even when the target is met |
| `"eq"`  | equal to `target`                   | squared relative deviation                                                                                             |


`weight` sets the relative importance. A target of zero or below has no
relative scale and needs `Spec(..., scale=...)` in the unit of the spec.

**Corners.** `Corner(name, table, nmos_name, pmos_name, conditions)`. The first
corner is the reference at which the design is sized.

**`run()` arguments.**

| Argument         | Default     | Meaning                                                             |
| ---              | ---         | ---                                                                 |
| `output_module`  | required    | path of the generated netlist                                       |
| `maxiter`        | `300`       | CMA-ES iterations                                                   |
| `n_restarts`     | `1`         | independent CMA-ES runs; raise it if results vary between seeds     |
| `seed`           | `1`         | random seed                                                         |
| `ref_index`      | `0`         | index of the reference corner                                       |
| `workers`        | `1`         | parallel processes, see [Parallel evaluation](#parallel-evaluation) |
| `netlist_format` | `"spectre"` | netlist format                                                      |


### Step 11: Run and read the results

```bash
python config.py
```

The run prints the CMA-ES progress, then the sized transistors (from the
reference corner):

```
Transistor Details:
  M1a:
    Length:  253.7nm
    Width:   5.473µm
    Area:    1.388 µm²
    Current: 5.561µA
    GmID:    19.98
  ...
```

followed by the performance at every corner:

```
Corner Results:
  Spec         | tt        | ff        | ss        | Target  | Binding
  -------------+-----------+-----------+-----------+---------+--------
  GBW          | 10.29M    | 10.47M    | 9.888M    | ≥10M    | ss
  AC Gain (dB) | 56.57     | 55.97     | 56.85     | ≥50     | ff
  PM           | 73.28     | 74.08     | 72.56     | ≥70     | ss
  DC CMR (dB)  | 57.84     | 57.42     | 57.88     | ≥50     | ff
  VOUT_Error   | 6.235µ    | 318.8µ    | 383.8µ    | ≤20m    | ss
  Area         | 40.87 µm² | 40.87 µm² | 40.87 µm² | ≤50 µm² | tt
  Itotal       | 49.09µ    | 50.26µ    | 47.85µ    | ≤50µ    | ff
  VBP          | 740.6m    | 773m      | 704m      |         |
  VOUT_MAX     | 1.044     | 1.046     | 1.042     |         |
  ...
```

- **Binding** names the corner with the worst value of each spec. That corner
  limits the design.
- **Area** is identical at every corner because the geometry is frozen.
  **Itotal** and **VBP** vary, because the bias current is conserved and the
  gate voltage adjusts to it.

The netlist for the reference corner is written to `output_module`:

```
parameters vcm=400m vout_dc=600m
subckt amp (VINN VINP VOUT vdd vss)
M1a (n01 VINN n02 vdd) pmos l=253.7n w=5.473u
M1b (n03 VINP n02 vdd) pmos l=253.7n w=5.473u
M2a (n01 n01 vss vss) nmos l=2u w=1.832u
M2b (n03 n01 vss vss) nmos l=2u w=1.832u
M3 (n02 vbp vdd vdd) pmos l=422.2n w=3.442u
M4 (VOUT n03 vss vss) nmos l=1.913u w=11.97u
M5 (VOUT vbp vdd vdd) pmos l=422.2n w=11.75u
Rz (n04 n03) resistor r=9.995k
CC (n04 VOUT) capacitor c=1.363p
Mvbp (vbp vbp vdd vdd) pmos l=422.2n w=3.442u
Ivbp (vbp vss) isource dc=11.29u
ends amp
```

The last two lines are the bias replica generated from
`VSource("VBP", mirror="M3")` and `mirror_currents()`. The subcircuit can be
placed directly in a testbench. Simulating it at every corner is the final
check of the model: the simulated operating point should match the reported
one.

## Checklist for a new circuit

1. Draw the schematic and name every node.
2. Transcribe the topology into `MOSFETS`, `PASSIVES`, and `VSOURCES`.
3. Choose the knobs: $\mathrm{gm}/I_D$ and length per matched group, one current per
   independent branch, ratios for mirrors, passive values.
4. Write every node voltage and every $V_{DS}$ on paper. Mark the
   self-referential ones as unknowns.
5. Write `solve_point()`: conditions, lookups in dependency order, currents,
   widths, lengths, $\mathrm{gm}/I_D$, small-signal parameters, `State`.
6. Write one residual per unknown.
7. Write `specs()`, including a centering error if an output voltage is imposed
   as a condition.
8. If a bias is generated by a mirror or derived from a knob, add the
   multicorner conservation hooks.
9. Add the netlist hooks for passives and biases.
10. Write the configuration file, run it, and simulate the netlist.


## Parallel evaluation

With `workers=4`, the candidates of each CMA-ES generation are evaluated by
four processes; `workers=1` (the default) evaluates them serially. The
processes share the lookup tables in memory, and the SLSQP polish runs
serially afterwards.

The processes are started with `spawn` on every platform. Therefore:

- keep the call to `run()` inside `if __name__ == "__main__":`;
- define the circuit class in an importable module (such as `design.py`), not
  in a notebook. For a class defined in a notebook, use `workers=1`.
