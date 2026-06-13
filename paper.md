---
title: 'Angelica: An open-source Python platform for pipe network simulation'
tags:
  - Python
  - pipe networks
  - hydraulics
  - oil and gas
  - flow simulation
  - pressure correction
authors:
  - name: Simon Rodriguez
    orcid: 0000-0001-6924-2428
    affiliation: 1
affiliations:
  - name: University College Dublin, Ireland
    index: 1
date: 12 June 2026
bibliography: paper.bib
---

# Summary

Angelica is an open-source Python platform for steady, incompressible, isothermal
simulation of pipe networks. It solves for nodal pressures and component flow rates
in networks of pipes, fittings, and pumps using a pressure-correction algorithm
(modified SIMPLE) with Colebrook-White [@colebrook1939], Hazen-Williams
[@williams1920], or Hagen-Poiseuille closures. A native graphical network editor
allows users to build, run, and inspect networks interactively without writing code.
Fluid properties for water and general liquids are entered directly; crude oil
properties — density and viscosity — are computed automatically from API gravity and
operating temperature via the Beggs & Robinson dead-oil correlation [@beggs1975].

# Statement of Need

Pipe network simulation sits between two ecosystems that have historically not
communicated: water distribution tools and oil & gas pipeline simulators.

On the water side, EPANET [@rossman2000] is the established open-source standard
for municipal networks. It is thoroughly validated but written in C; its Python
interface, WNTR [@klise2017], exposes EPANET's feature set through a scripting
API with no graphical editor. Neither tool provides oil & gas fluid property
support out of the box.

On the oil & gas side, commercial packages such as PIPESIM and PIPEPHASE handle
multiphase and compositional flows but are proprietary and expensive, creating a
barrier for research, education, and small-scale engineering projects.

Angelica addresses this gap as a Python-native platform with an accessible
graphical interface, a validated steady-state incompressible solver, and initial
oil & gas fluid property support. Its architecture is designed as a foundation for
future extensions toward non-isothermal, compressible, multiphase, and
compositional simulation — the capability set needed for full-scale oil & gas
pipeline engineering.

# Implementation

**Solver.** Angelica implements a two-phase pressure-correction procedure inspired
by the SIMPLE algorithm [@patankar1980]. A laminar seed (Hagen-Poiseuille closure)
initialises the pressure and flow fields; a turbulent correction phase then iterates
to convergence using the user-selected friction-factor closure. Two criteria must be
satisfied simultaneously for convergence: the maximum absolute pressure correction
must be below $10^{-3}$ Pa, and the maximum nodal mass-flow imbalance must be below
$10^{-3}$ relative to the total inflow. The linear system assembled at each
iteration is solved with a direct sparse solver via SciPy.

**Pressure-drop closures.** Three closures are available for pipes:

- *Colebrook-White* [@colebrook1939]: iterative solution of the implicit friction-factor equation for turbulent flow ($Re \geq 2300$); the laminar branch $f = 64/Re$ is applied automatically when $Re < 2300$, ensuring correct behaviour across the full Reynolds number range.
- *Hazen-Williams* [@williams1920]: empirical head-loss formula widely used in water distribution practice.
- *Hagen-Poiseuille*: exact analytical closure for fully developed laminar viscous flow, used as the initialisation seed for every simulation.

Fittings are modelled as local losses through user-defined or library K coefficients.
Pumps are represented by either a single-point head-flow model compatible with
EPANET conventions or a piecewise-linear multi-point characteristic curve.

**Fluid properties.** The `SingleComponentFluid` model accepts user-specified density
and viscosity for any single-phase liquid. The `dead_oil` module computes density
from API gravity via the standard specific-gravity conversion

$$\rho = \frac{141.5}{API + 131.5} \times 999.064 \; \text{kg m}^{-3}$$

and dynamic viscosity from the Beggs & Robinson dead-oil correlation [@beggs1975]

$$\log_{10}\!\left(\log_{10}(\mu_{\text{dead}} + 1)\right) = (3.0324 - 0.02023 \cdot API) - 1.163 \log_{10}(T_{\!\text{°F}})$$

where temperature appears only as the reference point for viscosity evaluation;
the simulation itself remains isothermal.

**Interface.** A Tkinter graphical editor provides drag-and-drop network
construction, real-time convergence monitoring, and tabular result export.
A self-contained Windows installer (no Python installation required) is
distributed alongside the pip-installable package, lowering the barrier to
adoption for practising engineers. Networks are serialised as JSON files that
can be loaded programmatically through the Python API.

# Validation

Eleven tutorial cases are distributed with Angelica; four serve as quantitative
benchmarks and one demonstrates the crude oil workflow:

- **Hagen-Poiseuille** (Tutorial 10): three parallel smooth pipes in the laminar
  regime ($Re < 180$) are solved and compared against the analytical Poiseuille
  solution $Q = \pi D^4 \Delta P / (128 \mu L)$. The maximum relative error across
  all pipes is below 0.001%.
- **Three-reservoir junction** (Tutorial 9): a classic three-reservoir problem
  solved with Colebrook-White friction. The computed junction head matches the
  published reference solution to within 0.01 m.
- **Hanoi network** (Tutorial 7): the 32-node, 34-pipe Hanoi benchmark network
  [@fujiwara1990]; flows and pressures are compared against EPANET [@rossman2000]
  reference output.
- **EPANET pump network** (Tutorial 8): a pump-and-pipe network validated against
  EPANET [@rossman2000] results.
- **Crude oil pipeline** (Tutorial 11): a branched gathering pipeline carrying
  32°API crude at 65°C. Fluid density (864.6 kg m$^{-3}$) and viscosity (4.25 cP)
  are computed automatically via the Beggs & Robinson correlation [@beggs1975]; all
  three pipes operate in the fully turbulent regime ($Re = 19\,000$–$31\,000$).

# Acknowledgements

The author thanks the open-source hydraulics community, whose published benchmarks
and reference solutions made rigorous validation possible.

# References
