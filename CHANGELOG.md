# Changelog

## [1.5.1] — 2026-08-17

### Tutorials
- Fix tutorial 02 GUI file: nodes 2 and 3 were at the same x-position,
  making the loop look like two horizontal parallel lines.  Repositioned
  to form a clear diamond (◇) shape.  Also corrected pipe lengths to
  match `run.py` (lower path: 8 km + 3 km, not 5.5 + 5.5 km).
- Add tutorial 03 — diamond loop with **flow outlet** BC: same 4-pipe
  diamond network as tutorial 02 but the outlet specifies mass-flow
  demand (100 kg/s) instead of pressure.  The solver determines the
  outlet pressure (5.356 MPa) and the flow split (63 / 37 %).
- Add tutorial 04 — looped network with **two flow delivery points**:
  5-node, 5-pipe network combining a diamond loop with a branch to a
  satellite separator.  Both outlets use flow BCs (80 kg/s + 20 kg/s).
  The loop routes ≈ 51 % via the upper path and ≈ 49 % via the lower
  path; all nodes remain undersaturated (P > 5.63 MPa bubble point).

## [1.5.0] — 2026-08-17

### New: Black-oil GUI integration

- Physics selector: new "Black-oil (3-phase)" option alongside Compressible
  and Incompressible modes.  Selecting it disables the energy radio buttons
  (black-oil always solves the energy equation).
- Material dialog: new "Black-oil" definition mode with four fields —
  API gravity, gas gravity, GOR sc (m³/m³), and WOR sc (m³/m³).
  No density or viscosity input required; PVT correlations compute them
  from first principles at run time.
- Node properties panel: thermal boundary condition section is now shown
  for black-oil cases (same as non-isothermal / compressible) — at least
  one "Fixed temperature" node is required.
- Link properties panel: heat-transfer fields (U and T_amb) and the
  "Heat Source" component button are now available in black-oil mode.
- Two GUI tutorial files added:
  - `steady_black_oil/01_three_phase_pipeline/three_phase_pipeline.gui.json`
    — 10 km, 0.20 m pipe; 8 MPa inlet (undersaturated) → 2 MPa outlet
    (two-phase below bubble point ~5.63 MPa).
  - `steady_black_oil/02_looped_gathering_network/looped_gathering_network.gui.json`
    — 4-node, 4-pipe loop; two parallel paths from 8 MPa to 2 MPa;
    demonstrates split-flow with phase-state differences between paths.

## [1.4.2] — 2026-08-17

### Fixes
- All three thermal solvers (non-isothermal, compressible, black-oil) now raise a
  clear `ValueError` at the start of `solve()` if no `fixed_temperature` thermal
  boundary is provided. Without a Dirichlet node the energy system is singular
  (pure Neumann problem) and previously produced silently incorrect temperatures.

## [1.4.1] — 2026-08-17

### Fixes
- Fix compressible solver: inlet temperature boundary condition was ignored during
  solving. The solver seeded all nodes with the inlet temperature but never marked
  the inlet node as a Dirichlet boundary (`is_thermal_inlet = True`), so the energy
  solver was free to change it every iteration. The prescribed temperature is now
  locked correctly, matching the behaviour of the non-isothermal solver.

## [1.4.0] — 2026-08-16

### New: Black-oil three-phase solver

- `BlackOilFluid` — homogeneous three-phase fluid model (stock-tank oil +
  dissolved/free gas + water) implementing the `FluidModel` interface.
  Works with all existing network components without changes to the core
  hydraulic or energy solvers.
- `SteadyBlackOilSolver` — steady-state solver with outer PVT iteration loop.
  Converges when the mixture-density field and temperature field both stop
  changing between iterations.
- PVT correlations:
  - Bubble point, Rs, Bo: Standing (1947)
  - Gas compressibility z: Hall-Yarborough (1974) + Sutton (1985) pseudo-crits
  - Gas viscosity: Lee, Gonzalez & Eakin (1966)
  - Live-oil viscosity: Beggs & Robinson (1975)
  - Water FVF: McCain (1990), simplified
- 31 new unit tests covering all correlations and phase-holdup logic.
- Tutorial `steady_black_oil/01_three_phase_pipeline`: 32°API crude, GOR 25 m³/m³,
  WOR 0.5, single 10 km pipe — demonstrates undersaturated inlet transitioning
  to two-phase at the outlet when pressure drops below the bubble point.

## [1.3.31] — 2026-08-16

### Solver improvements
- Fix energy solver accuracy: use the exact analytical NTU formula
  (`T_exit = T_amb + (T_in − T_amb)·exp(−NTU)`) for pipes with
  `n_thermal_segments = 1`, replacing the single-cell FV approximation.
- Change `velocity_loop_method` default from `"fixed_point"` to `"secant"`.
  The secant method converges robustly from any initial guess, including
  all-pressure-boundary cases where the outer loops exit after one iteration
  and the fixed-point iteration had insufficient evaluations.

### GUI
- Fix re-run simulation bug (results panel was not reset between runs).

### Packaging
- `sv-ttk` moved to optional `[gui]` extra — headless API use no longer
  requires the Tkinter theme package.
- `openpyxl` moved to optional `[excel]` extra — solver-only installs no
  longer pull in a spreadsheet library.
- Hard dependencies reduced to `numpy` and `scipy`.
- Added `authors`, `keywords`, trove `classifiers`, and `[project.urls]`
  to `pyproject.toml`.
- Added `py.typed` marker; package is now typed per PEP 561.
- Added `pytest-cov`, `ruff`, `mypy` to the `[dev]` optional extra.

### Type annotations
- `BaseSolver.solve(case: NetworkCase, ...) -> SolveResult` annotated.
- Annotation propagated to all three concrete solver `solve()` methods.

### Code quality
- Extracted `_update_component_temperatures` from both outer solvers into
  `BaseSolver`.
- Extracted `_initial_temperature` from both outer solvers into `BaseSolver`
  (the non-isothermal version had a spurious unused `_network_state` parameter).
- Deleted legacy shims `angelica/solver.py` and `angelica/reporting.py`.
- `export_solve_result_workbook` now raises a clear `ImportError` with install
  instructions when `openpyxl` is not available.

### Tests
- Added three quantitative compressible benchmarks:
  - P²-law mass flow (single pipe vs. analytical P₁²−P₂² + Colebrook-White).
  - Series intermediate pressure (two isothermal pipes, P_mid = √((P₁²+P₂²)/2)).
  - NTU heat loss (compressible pipe vs. analytical exit temperature).
- Added flow-value assertions to six previously convergence-only tutorial tests.
- Fixed `test_reporting.py` to skip gracefully when `openpyxl` is not installed
  (`pytest.importorskip`).

### Tutorials
- Added `run.py` and `README.md` to compressible tutorial 02
  (flow BC cross-validation — mixed pressure/flow boundaries).
- Added `run.py` to compressible tutorial 03 (looped gas pipeline with heat
  loss — 7-pipe methane network, 700→500 kPa).
- Updated tutorial 03 compressible `README.md` with actual solver output.
- Added `README.md` to non-isothermal tutorial 05 (32°API crude oil pipeline,
  temperature-dependent viscosity, 44× cold/hot viscosity ratio).
- Redesigned non-isothermal looped network tutorial for physical clarity.

### Documentation
- `README.md`: compressible solver marked as live (✓) in the roadmap graphic.
- `README.md`: added compressible row to features table; tutorials section
  expanded to 19 cases across three solver families.
- `tutorials/README.md`: added tutorial 05 non-isothermal and all three
  compressible tutorials.


## [1.3.30] — 2026-08-15

### New
- Compressible tutorial 02: flow BC cross-validation.
- Compressible tutorial 03: looped gas pipeline with heat loss (7-pipe network
  with feeder, upper/lower paths, diagonal cross-pipe, and collector).

### Fixes
- Fix listbox deselection closing the component properties panel.


## [1.3.29] — 2026-08-14

### New
- Pipe thermal segmentation: `num_segments` field exposed in the GUI, allowing
  per-pipe control over the number of finite-volume energy cells.


## [1.3.28] — 2026-08-14

### Fixes
- Compressible solver now requires at least one fixed-temperature thermal
  boundary; raises a clear error if none is provided.


## [1.3.27] — 2026-08-13

### Fixes
- Fix compressible tutorial 01 inlet thermal boundary condition.
