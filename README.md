<div align="center">

# Angelica

**The open-source platform for pipeline network simulation.**

[![License: MIT](https://img.shields.io/badge/License-MIT-0c6c6b.svg)](LICENSE)
[![Python](https://img.shields.io/badge/Python-3.8%2B-3776ab)](https://www.python.org/)
[![Website](https://img.shields.io/badge/Website-scimonCFD.github.io-d98f3b)](https://scimonCFD.github.io/Angelica_website/)

</div>

Angelica solves **steady, incompressible, isothermal pipe networks** — computing nodal pressures and component flow rates for systems with pipes, fittings, and pumps. A graphical interface means there is nothing to script: define sources and sinks, connect them through junctions, assign components, and run.

## Features

| Component | Description |
|-----------|-------------|
| **Pipes** | Darcy-Weisbach with Colebrook-White or Hazen-Williams closure. Arbitrary length, diameter, and wall roughness. |
| **Fittings** | Local-loss accessories from a built-in library or a custom K coefficient. |
| **Pumps** | Single-point EPANET model and piecewise-linear multi-point curves. |
| **Boundaries** | Pressure and mass-flow boundaries on any node. Mixed types on the same network. |
| **GUI** | Graphical network editor with convergence monitor and spreadsheet report export. |
| **Validation** | Results benchmarked against published EPANET reference solutions. |

## Installation

```bash
git clone https://github.com/ScimonCFD/Angelica.git
cd Angelica
pip install -e .
```

Requires Python 3.8+ and the dependencies listed in [`pyproject.toml`](pyproject.toml).

## Quick Start

Launch the graphical interface:

```bash
python -m angelica.gui.app
```

**Workflow:**
1. Add **source** and **sink** nodes
2. Connect them through **junction** nodes as needed
3. Assign **pipes, fittings, or pumps** to each link
4. Set fluid properties and solver settings
5. **Run** — inspect convergence, review nodal and link results, export a report

## Tutorials and Benchmarks

Eight cases are included under [`tutorials/steady_isothermal_incompressible/`](tutorials/steady_isothermal_incompressible/):

| # | Case | Type |
|---|------|------|
| 01 | Pipe-only base case | Tutorial |
| 02 | Network with fittings | Tutorial |
| 03 | Fittings and elevation changes | Tutorial |
| 04 | Inlet mass-flow boundary | Tutorial |
| 05 | Outlet mass-flow boundary | Tutorial |
| 06 | Combined inlet and outlet mass-flow boundaries | Tutorial |
| 07 | Hanoi network — 32 nodes, 34 pipes | Benchmark (EPANET reference) |
| 08 | EPANET pump network | Benchmark (EPANET reference) |

## Roadmap

The incompressible solver is the first stage of a platform designed to grow:

```
Steady isothermal  →  Non-isothermal  →  Compressible  →  Multiphase  →  Multicomponent
                                                                               ↑
                                                               Long-term: Transient
```

## Repository Layout

```
src/angelica/
├── core/        # Network topology, components, state, settings, results
├── properties/  # Fluid-property models
├── closures/    # Pressure-drop and device models
├── numerics/    # Assembly, convergence, and linear algebra
├── solvers/     # Solver implementations
├── io/          # Reporting helpers
└── gui/         # Graphical network editor
```

## Documentation

- **[Methods page](https://scimonCFD.github.io/Angelica_website/methods.html)** — pressure-correction algorithm and element laws
- **[`docs/WORKFLOW_AND_EQUATIONS.md`](docs/WORKFLOW_AND_EQUATIONS.md)** — governing equations and solver notes

## Contact

Angelica is developed by [Simon Rodriguez](https://www.linkedin.com/in/simonrodriguezl), Venezuelan engineer and researcher based in Dublin, Ireland.

Questions, suggestions, or bug reports — [open an issue](https://github.com/ScimonCFD/Angelica/issues) or reach out on [LinkedIn](https://www.linkedin.com/in/simonrodriguezl).

## License

MIT — see [LICENSE](LICENSE).
