<div align="center">

# Angelica

**The open-source platform for pipeline network simulation.**

[![License: MIT](https://img.shields.io/badge/License-MIT-0c6c6b.svg)](LICENSE)
[![Python](https://img.shields.io/badge/Python-3.8%2B-3776ab)](https://www.python.org/)
[![Website](https://img.shields.io/badge/Website-scimonCFD.github.io-d98f3b)](https://scimonCFD.github.io/Angelica_website/)

</div>

Angelica solves **steady, incompressible pipe networks** — computing nodal pressures, flow rates, and (in non-isothermal mode) temperatures for systems with pipes, fittings, pumps, and heat sources. A graphical interface means there is nothing to script: define sources and sinks, connect them through junctions, assign components, and run.

## Features

| Component | Description |
|-----------|-------------|
| **Pipes** | Darcy-Weisbach with Colebrook-White or Hazen-Williams closure. Arbitrary length, diameter, and wall roughness. |
| **Fittings** | Local-loss accessories from a built-in library or a custom K coefficient. |
| **Pumps** | Single-point EPANET model and piecewise-linear multi-point curves. |
| **Heat sources** | Inline heaters and coolers with fixed power (W) and optional pressure drop. Non-isothermal mode only. |
| **Boundaries** | Pressure and mass-flow boundaries on any node. Mixed types on the same network. |
| **Non-isothermal** | Outer temperature loop with NTU-based pipe heat loss, source-term linearisation (Moukalled), and temperature-dependent fluid properties. |
| **GUI** | Graphical network editor with convergence monitor and spreadsheet report export. |
| **Validation** | Results benchmarked against published EPANET reference solutions and NTU analytical solutions. |

## Installation

### Windows — one-click installer

Download **AngelicaSetup.exe** from the [latest release](https://github.com/ScimonCFD/Angelica/releases/latest) and run it. No Python required. The installer places Angelica in your user folder (no admin rights needed) and optionally creates a desktop shortcut.

### Linux / macOS — pip from GitHub

Requires Python 3.8 or later.

```bash
pip install git+https://github.com/ScimonCFD/Angelica.git
```

Then launch the GUI:

```bash
angelica-gui
```

### Development install

```bash
git clone https://github.com/ScimonCFD/Angelica.git
cd Angelica
pip install -e .
angelica-gui
```

## Quick Start

Launch the graphical interface:

```bash
angelica-gui
```

**Workflow:**
1. Add **source** and **sink** nodes
2. Connect them through **junction** nodes as needed
3. Assign **pipes, fittings, or pumps** to each link
4. Set fluid properties and solver settings
5. **Run** — inspect convergence, review nodal and link results, export a report

## Tutorials and Benchmarks

Fifteen cases are included across two solver folders.

**Isothermal** — [`tutorials/steady_isothermal_incompressible/`](tutorials/steady_isothermal_incompressible/)

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
| 09 | Three-reservoir junction | Benchmark (analytical Colebrook-White) |
| 10 | Laminar Poiseuille parallel pipes | Benchmark (exact closed-form) |
| 11 | Crude oil gathering pipeline — 32°API, 65°C | Demo (Beggs & Robinson properties) |

**Non-isothermal** — [`tutorials/steady_non_isothermal_incompressible/`](tutorials/steady_non_isothermal_incompressible/)

| # | Case | Type |
|---|------|------|
| 12 | Single pipe with heat loss to ambient | Benchmark (NTU analytical) |
| 13 | Branched district-heating network with two thermal loads | Tutorial |
| 14 | Looped network with ambient heat loss | Tutorial |
| 15 | Inline electric heater — energy balance verification | Benchmark (ΔT = Q/ṁcₚ) |

## Roadmap

Two solver modes are live; the platform is designed to grow:

```
Steady isothermal ✓  →  Non-isothermal ✓  →  Compressible  →  Multiphase  →  Multicomponent
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
