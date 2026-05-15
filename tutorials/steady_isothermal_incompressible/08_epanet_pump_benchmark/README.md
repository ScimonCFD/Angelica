# EPANET Pump Benchmark

This folder contains a pump-bearing validation case derived from the official
EPANET 2.2 tutorial network.

It is useful for `netSim` because it exercises two things at once:

- looped-network hydraulics with elevation changes
- a pump represented through a head-versus-flow curve

## Reference

This case is adapted from the official EPANET 2.2 documentation:

- `2. Quick Start Tutorial`, which defines the tutorial network and the pump
  design point
- `Appendix C`, which publishes the exact `.inp` definition and the reference
  report values at `0:00`

Official source:

- https://usepa.github.io/EPANET2.2/2_quickstart.html
- https://usepa.github.io/EPANET2.2/back_matter.html

## What the Original Case Contains

The original EPANET tutorial network contains:

- one reservoir
- one pump between the reservoir and the main network
- a two-loop pipe network
- one elevated tank
- four demand nodes downstream of the pumped system

The EPANET tutorial input used here uses:

- `Hazen-Williams`
- `GPM` flow units
- a one-point pump curve:
  - `Q = 1000 gpm`
  - `Head = 200 ft`

EPANET expands that one-point curve internally into a complete pump curve.

## How It Is Represented in `netSim`

The current `netSim` node model keeps junctions conservative, so the four
published nodal demands are represented as side branches to dedicated `sink`
nodes, just as in the validated Hanoi adaptation.

The benchmark comparison is performed on the original trunk network only:

- trunk-node pressures
- trunk-node total heads
- trunk-link flows

The side-branch demand sinks are only a modeling device used to fit the EPANET
tutorial into the current `netSim` node model. They are not themselves part of
the published benchmark tables.

## Pump Modeling

In this case, the pump is entered as a `Q-Head` table and solved in `netSim`
as another pressure-changing component:

- the current pressure field implies a pump differential head
- the pump curve is inverted locally to obtain `Q`
- the local derivative of the curve supplies the pressure-correction coupling

For the official EPANET one-point curve, `netSim` uses the same EPANET
completion rule:

- shutoff head at zero flow = `133%` of the design head
- zero head at twice the design flow

## Expected Results

The published EPANET report at `0:00` provides:

- node heads and pressures
- link flows

The script
[benchmark_tables.py](/home/simon/Documents/netSim/tesisIca/ClasesTesis/netSim/netSim/tutorials/steady_isothermal_incompressible/08_epanet_pump_benchmark/benchmark_tables.py)
prints:

- expected node pressures vs `netSim`
- expected node total heads vs `netSim`
- expected trunk-link flows vs `netSim`

## Useful Commands

Run the tutorial case:

```bash
python3 tutorials/steady_isothermal_incompressible/08_epanet_pump_benchmark/run.py
```

Print the benchmark comparison tables:

```bash
python3 tutorials/steady_isothermal_incompressible/08_epanet_pump_benchmark/benchmark_tables.py
```

Open the case in the GUI:

```bash
PYTHONPATH=src python3 -m netSim.gui.app
```

then use:

- `File -> Open`
- open
  [epanet_pump_benchmark.gui.json](/home/simon/Documents/netSim/tesisIca/ClasesTesis/netSim/netSim/tutorials/steady_isothermal_incompressible/08_epanet_pump_benchmark/epanet_pump_benchmark.gui.json)
- press `Run`
