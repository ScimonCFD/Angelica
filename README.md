# netSim

**netSim** is an open-source hydraulic network simulator for **steady-state, incompressible, single-phase flow**.

Its current release focuses on a narrow but useful core: solving network pressures and flow rates in pipes, fittings, and pumps, with a GUI-first workflow and a codebase designed to grow toward more general pipeline simulation.

The long-term direction is broader than water networks alone. The project is being built as a modular foundation for future extensions toward **oil, gas, natural gas, and black-oil systems**, while keeping the current solver lightweight, inspectable, and easy to extend.

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

## Why This Project Exists

Many established tools are strong inside their original domain, but harder to extend when you want to move beyond a classical water-network workflow.

`netSim` is meant to be different:

- transparent enough to study and modify
- structured enough to validate and maintain
- flexible enough to evolve toward richer physics over time

That makes it useful both as:

- a small but real hydraulic network solver
- and a foundation for future research and development

## Validation

The current solver is not just running toy cases. It has already been checked against recognized benchmark-style references, including:

- the **Hanoi** EPANET/Darcy benchmark network
- an **EPANET pump tutorial benchmark** with published reference node pressures, heads, and link flows

These cases live in the tutorial suite and can be run directly from the repository.

## Quick Start

From the repository root:

```bash
pip install -e .
```

Then define and solve a minimal case:

```python
from netSim.core.case import NetworkCase, PressureBoundary
from netSim.core.components import Pipe
from netSim.properties.single_component import SingleComponentFluid
from netSim.solvers.steady_isothermal_incompressible import SteadyIsothermalIncompressibleSolver

fluid = SingleComponentFluid(
    density_kg_per_m3=998.25,
    viscosity_pa_s=0.001,
)

pipe = Pipe(
    start_node=1,
    end_node=2,
    length_m=100.0,
    diameter_m=0.1,
    absolute_roughness_m=0.000046,
)

case = NetworkCase(
    name="example",
    fluid_model=fluid,
    components=[pipe],
    pressure_inlets=[PressureBoundary(node_id=1, pressure_pa=200_000.0)],
    pressure_outlets=[PressureBoundary(node_id=2, pressure_pa=100_000.0)],
)

solver = SteadyIsothermalIncompressibleSolver()
result = solver.solve(case)
print(result)
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

Validated examples are available under
`tutorials/steady_isothermal_incompressible/`.

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

This version is the **foundation release**, not the final target.

Planned next steps include:

- richer liquid-property modes, including oil-oriented workflows
- non-isothermal solving through an energy equation
- compressible single-phase flow
- more active devices and operating logic
- black-oil and broader hydrocarbon-network support

The goal is for those future capabilities to build on the same network core rather than replace it.

## Technical Notes

The solver workflow and governing-equation notes are documented in:

[`docs/WORKFLOW_AND_EQUATIONS.md`](/home/simon/Documents/netSim/tesisIca/ClasesTesis/netSim/netSim/docs/WORKFLOW_AND_EQUATIONS.md)
