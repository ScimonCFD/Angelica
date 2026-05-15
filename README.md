# netSim

**netSim** is an open-source hydraulic network simulator for **steady-state, incompressible, single-phase flow**.

It computes nodal pressures and component flow rates in pipe networks with support for pipes, fittings, pumps, elevation changes, and pressure or mass-flow boundary conditions.

## Current Scope

Today, `netSim` supports:

- steady-state, incompressible, isothermal flow
- laminar and turbulent pipe flow
- pipe pressure-drop models:
  - Darcy-Weisbach with Colebrook-White
  - Hazen-Williams
- fittings via local-loss coefficients
  - including a library of named preset components
  - plus manual user-defined `K`
- pumps defined by `Q-Head` tables
- elevation changes through a gravitational pressure term
- pressure boundaries and mass-flow boundaries
- a graphical editor for building and running cases
- export of converged results to spreadsheet report files

The numerical core uses a **segregated pressure-correction method** with adaptive laminar initialisation and explicit pressure relaxation.

## Validation

The current solver has been checked against benchmark-style reference cases, including:

- the **Hanoi** EPANET/Darcy benchmark network
- an **EPANET pump tutorial benchmark** with published reference node pressures, heads, and link flows

These cases live in the tutorial suite and can be run directly from the repository.

## Quick Start

From the repository root:

```bash
pip install -e .
```

Launch the GUI:

```bash
python -m netSim.gui.app
```

Then:

1. place `source`, `sink`, and `junction` nodes
2. connect them with links
3. add `pipe`, `fitting`, or `pump` components to each link
4. define the material
5. choose the pipe pressure-drop model and numerics
6. press `Run`

To start from a validated example, open one of the tutorial cases from:

- [tutorials/steady_isothermal_incompressible](/home/simon/Documents/netSim/tesisIca/ClasesTesis/netSim/netSim/tutorials/steady_isothermal_incompressible)

For example:

```bash
python tutorials/steady_isothermal_incompressible/07_hanoi_epanet_darcy_benchmark/run.py
python tutorials/steady_isothermal_incompressible/08_epanet_pump_benchmark/run.py
```

## GUI

A graphical network editor is included.

Launch it with:

```bash
python -m netSim.gui.app
```

The GUI currently supports:

- visual network construction
- source, sink, and junction placement
- pipe, fitting, and pump definition
- solver and correlation configuration
- convergence inspection
- report export for converged runs

## Tutorials

Validated examples are available under `tutorials/steady_isothermal_incompressible/`.

| # | Case |
|---|------|
| 01 | Pipe-only base case |
| 02 | Network with fittings |
| 03 | Network with fittings and elevation changes |
| 04 | Inlet mass-flow boundary |
| 05 | Outlet mass-flow boundary |
| 06 | Combined inlet and outlet mass-flow boundaries |
| 07 | Hanoi benchmark (EPANET reference) |
| 08 | EPANET pump benchmark |

Run a case, for example:

```bash
python tutorials/steady_isothermal_incompressible/07_hanoi_epanet_darcy_benchmark/run.py
python tutorials/steady_isothermal_incompressible/08_epanet_pump_benchmark/run.py
```

## Repository Layout

```text
src/netSim/
├── core/        # Network topology, components, state, settings, results
├── properties/  # Fluid-property models
├── closures/    # Pressure-drop and device models
├── numerics/    # Assembly, convergence, and linear algebra helpers
├── solvers/     # Solver implementations
├── cases/       # Reusable case definitions
├── io/          # Reporting helpers
└── gui/         # Graphical network editor
```

## Roadmap

Planned next steps include:

- richer liquid-property modes, including oil-oriented workflows
- non-isothermal solving through an energy equation
- compressible single-phase flow
- more active devices and operating logic
- black-oil and broader hydrocarbon-network support

These future capabilities are intended to build on the same network core.

## Technical Notes

The solver workflow and governing-equation notes are documented in:

[`docs/WORKFLOW_AND_EQUATIONS.md`](/home/simon/Documents/netSim/tesisIca/ClasesTesis/netSim/netSim/docs/WORKFLOW_AND_EQUATIONS.md)
