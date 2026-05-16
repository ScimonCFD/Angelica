# netSim

`netSim` is a simulator for steady-state pipe networks with pipes, fittings, and pumps.

It computes nodal pressures and component flow rates for incompressible, isothermal, single-phase systems.

## Quick Start

`netSim` comes with a graphical user interface.

Launch it from the repository root with:

```bash
python -m netSim.gui.app
```

The workflow follows the traditional style used in pipe-network simulators:

1. define sources and sinks
2. connect them through internal junctions as needed
3. assign pipes, fittings, and pumps to the links
4. define the material and numerical settings
5. run the case from the GUI

The intended entry point for most users is the GUI and the user manual.

## Tutorials

Example and benchmark cases are available under
[tutorials/steady_isothermal_incompressible](/home/simon/Documents/netSim/tesisIca/ClasesTesis/netSim/netSim/tutorials/steady_isothermal_incompressible).

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

These cases are useful both as examples and as validation references for the current solver.

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

## Long-Term Objective

Planned next steps include:

- steady, non-isothermal flow
- multiphase flow
- multicomponent flow
- broader hydrocarbon-network support

## Publications

To be defined.

## Contact

Simon Rodriguez

## Technical Notes

The solver workflow and governing-equation notes are documented in:

[`docs/WORKFLOW_AND_EQUATIONS.md`](/home/simon/Documents/netSim/tesisIca/ClasesTesis/netSim/netSim/docs/WORKFLOW_AND_EQUATIONS.md)
