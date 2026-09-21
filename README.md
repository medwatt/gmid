# Mosplot

**Mosplot** is a Python framework for analog circuit design based on the **$\mathrm{gm}/I_D$
methodology**. It characterizes MOS transistors through circuit simulation,
provides design charts derived from the resulting data, and performs automated
multicorner sizing of complete circuits.

![Examples of charts produced with Mosplot](docs/figures/overview.svg)

## Overview

The framework consists of three components.

### 1. Lookup table generation

A DC sweep over channel length and terminal voltages is performed with
**ngspice**, **HSPICE**, or **Spectre**. The operating-point parameters of the
device (drain current, transconductances, output conductance, and capacitances)
are stored in a single lookup table. One table is generated per technology and
process corner.

The simulator configuration, sweep definition, supported parameters, and
generation of corner tables are described in [Lookup table
generation](docs/lookup_tables.md).

### 2. Design charts

The conventional $\mathrm{gm}/I_D$ design flow is supported directly. Current
density, transit frequency, intrinsic gain, or any user-defined expression can
be plotted against $\mathrm{gm}/I_D$ for each channel length, and exact values
at any operating point are obtained by interpolation on the lookup table.

```python
nmos = Mosfet(lookup_table=lookup_table, mos="nmos", vbs=0.0, vds=0.6)
nmos.plot_by_expression(
    x_expression=nmos.gmid_expression,
    y_expression=nmos.transit_frequency_expression,
    y2_expression=nmos.gain_expression,
)
```

The available expressions and plotting functions are described in
[Design charts](docs/plotting.md).

### 3. Circuit optimization

A circuit is described once as a system of equations formulated in terms of
lookups into the pre-characterized table. The optimizer sizes the circuit
against the target specifications at all process corners simultaneously, reports
the performance at each corner, and generates a netlist for verification by
simulation.

Ready-made designs are included for single-ended and fully differential
amplifiers, current mirrors, and reference circuits. Using one of them in a new
technology requires only a configuration file with the lookup tables, operating
conditions, and target specifications. The list of designs, their use, and the
procedure for writing a new design are described in
[Circuit optimization](docs/optimization.md).

## Installation

Python 3.9 or later is required. Lookup table generation additionally requires
one of the supported simulators (ngspice, HSPICE, or Spectre).

```bash
pip install git+https://github.com/medwatt/gmid.git
```

For development, install from a local clone:

```bash
git clone https://github.com/medwatt/gmid.git
cd gmid
pip install -e .
```

## Citation

If this work is used in academic research, please cite:

```bibtex
@article{watfa2026residual,
  author  = {Watfa, Mohamed and Garcia-Ortiz, Alberto and Sassatelli, Gilles},
  title   = {A Residual-Equation Framework for {$\mathrm{gm}/I_D$} Analog Circuit Design
             with Automated Multicorner Sizing},
  journal = {IEEE Transactions on Circuits and Systems I: Regular Papers},
  year    = {2026},
  doi     = {10.1109/TCSI.2026.3736902}
}
```

## Acknowledgments

The HSPICE output parser is based on
[hspiceParser](https://github.com/HMC-ACE/hspiceParser).
