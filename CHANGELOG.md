# Changelog

## [1.6.44] — 2026-08-21

### Fix: Windows installer — "No module named pandas"

`pandas` was explicitly excluded from the PyInstaller bundle.  One or more
modules inside `thermo` (collected via `collect_all`) import pandas at runtime,
causing an immediate crash on Windows.  Removed `pandas` from the `excludes`
list so PyInstaller bundles it automatically alongside the other dependencies.

## [1.6.43] — 2026-08-21

### Change: Compositions moved to results report instead of canvas

Outlet and inlet compositions are no longer printed on the main canvas.
With mixtures of many components the canvas labels became cluttered and
unreadable.

Compositions now appear in the **results report** (CSV and Excel):

- **CSV** — a "Compositions" section at the end of the file with one row
  per pipe (mole fractions) and one row per node (mass-flow-weighted
  mixture at that junction).
- **Excel** — a dedicated **Compositions** sheet with the same data,
  formatted with a header row.

`SolveResult` carries `component_names` and `node_compositions` (both
computed by the compositional solver); the reporting functions use them
to write the new section/sheet automatically.  Non-compositional runs
produce no compositions section.

## [1.6.42] — 2026-08-21

### Feature: Stop button in Convergence window

A **Stop** button now appears in the top-right corner of the Convergence
Metrics window.  It is enabled while a simulation is running and disabled
otherwise.  Pressing it signals the solver thread to abort at the next
iteration callback, so heavy simulations can be cancelled without closing
the application.  The status bar shows "Simulation stopped by user."

### Feature: per-component composition tab for inlet nodes (compositional mode)

The source-node dialog in compositional mode now has two tabs:

- **Properties** — boundary type (pressure / flow), thermal BC (unchanged).
- **Composition** — one row per component as defined in the material:

  ```
  Component     Mole fraction
  ─────────────────────────────
  methane       [ 0.8 ]
  ethane        [ 0.2 ]
  propane       [     ]   ← blank → 0
  ```

  Blank fields are treated as 0.  The CSV `zs` string saved to the scene is
  built automatically from the individual entries when the user clicks Save.
  Old files with a pre-existing CSV value are pre-populated correctly.

## [1.6.41] — 2026-08-21

### Feature: outlet compositions shown on the canvas after a compositional solve

After running the compositional solver, each node in the network now shows
the local mole fractions alongside pressure, temperature, and flow rate:

- **Source nodes** — displays the inlet composition defined by the user
  (e.g. `CH₄=80.0%  C₂H₆=20.0%`) even before solving.
- **Sink and junction nodes** — displays the arriving/mixed composition
  computed by the solver (mass-flow-weighted average at junctions with
  multiple inlets).

**Implementation details:**

- `ComponentFlowResult` gains a new `zs: tuple[float, ...]` field (empty for
  non-compositional solvers, fully backward-compatible).
- `SteadyCompositionalSolver` populates `zs` from the converged per-pipe
  mole fractions.
- `_build_boundary_results` derives per-node compositions from the
  per-pipe `zs` values.
- `_node_summary_text` calls the new `_fmt_composition` helper to render
  compositions as `"CH₄=80.0%  C₂H₆=20.0%"`.
- Serialisation (`gui/io.py`) preserves `zs` in saved `.gui.json` files.

## [1.6.40] — 2026-08-21

### Fix: convergence plot markers clipped on the X axis too

v1.6.39 fixed vertical overflow of markers; this release fixes horizontal
overflow. When the first outer-pass boundary fell very close to `x = left`
(e.g. only one turbulent iteration in the first outer pass), the secondary-axis
diamond markers extended up to `size = 4 px` outside the left border.
Primary-series circle markers at the first or last iteration had the same
2 px overflow on the left/right borders.

Both marker types now clamp their pixel centres in both axes:
- Diamond: `px_m ∈ [left+4, right-4]`, `py_m ∈ [top+4, bottom-4]`
- Circle:  `xm  ∈ [left+2, right-2]`, `ym  ∈ [top+2, bottom-2]`

## [1.6.39] — 2026-08-21

### Fix: convergence plot markers no longer draw outside the plot area

Primary-series circle markers and secondary-series diamond markers could
extend a few pixels beyond the top or bottom boundary of the plot when a
data value landed exactly at the axis limits.

- Primary series: `y` pixel coordinate is now clamped to `[top, bottom]`
  before drawing both the polyline and the circle markers.
- Secondary series (right axis): diamond marker centres are clamped to
  `[top + size, bottom - size]` so the full diamond always stays inside
  the plot rectangle, regardless of where the value sits on the right axis.

## [1.6.38] — 2026-08-21

### Fix: Windows installer now bundles `thermo`, `chemicals`, and `fluids`

`angelica.spec` updated to use `collect_all("thermo")`, `collect_all("chemicals")`,
and `collect_all("fluids")`.  These packages use lazy imports and ship data
files (property databases) that PyInstaller cannot discover through static
analysis alone.  Without this fix the Windows `.exe` would raise
`ModuleNotFoundError: No module named 'thermo'` when opening a compositional case.

`build-windows.yml` already installs `thermo` automatically (it became a core
dependency in v1.6.37), so no workflow change is needed — only the spec changes.

## [1.6.37] — 2026-08-21

### Change: `thermo` is now a required dependency

`thermo>=0.5` moved from `[compositional]` optional extra to the core
`dependencies` list in `pyproject.toml`.  A plain `pip install angelica`
now installs `thermo` automatically — no extra flag needed.

The `[compositional]` extra is removed.  Users who had
`pip install "angelica[compositional]"` in their workflows can drop the
extra; `pip install angelica` is sufficient.

CI updated to install `.[dev]` (thermo comes in via core deps).
Skip markers removed from compositional tests — they run unconditionally.

## [1.6.36] — 2026-08-21

### Feature: GUI support for compositional physics mode

The desktop GUI now supports the **Compositional (EOS)** physics mode end-to-end:

**Physics → Case Type…**
- New "Compositional (EOS)" radio button in the Compressibility row.
- Energy row is disabled (compositional always solves the energy equation).

**Material → Define Material**
- Opens a dedicated compositional fluid dialog: comma-separated component
  names and default mole fractions (validated to sum to 1.0).
- Material summary panel shows component list and default zs.

**Node dialog (source nodes)**
- New "Inlet composition" section with a CSV mole-fractions field (zs).
- Required for every source node in compositional mode.
- Thermal boundary section is shown (same as non-isothermal/compressible).

**Pipe component properties**
- Heat-transfer coefficient and ambient temperature fields appear when the
  mode is compositional (same as other thermal modes).
- Heat Source component is available in the palette.

**Solver**
- `build_solver_from_scene` returns `SteadyCompositionalSolver` when
  `physics_mode == "compositional"`.
- `build_network_case_from_scene` constructs `CompositionalFluid` and
  `InletCompositionBC` objects from the scene data.

**Tutorial `.gui.json` files (NEW)**
- `tutorials/steady_compositional/01_single_pipe/01_single_pipe.gui.json`
- `tutorials/steady_compositional/02_gas_mixing_junction/02_gas_mixing_junction.gui.json`
- `tutorials/steady_compositional/03_looped_network/03_looped_network.gui.json`

All three open, build, and solve from the GUI without writing any Python.

## [1.6.35] — 2026-08-21

### Docs: README files for the three compositional tutorials

Added `README.md` to each compositional tutorial folder:

- `tutorials/steady_compositional/01_single_pipe/README.md`
- `tutorials/steady_compositional/02_gas_mixing_junction/README.md`
- `tutorials/steady_compositional/03_looped_network/README.md`

Each README documents the network topology, fluid definition, expected
output, and the key physics to observe.  The GUI does not yet support the
compositional physics mode — these tutorials are run as Python scripts.

## [1.6.34] — 2026-08-21

### Tutorial: compositional looped gas network

New tutorial `tutorials/steady_compositional/03_looped_network/run.py` and
case builder `cases/looped_gas_gathering.py`.

Two gas sources at different pressures and compositions feed a looped
five-pipe gathering network (Nodes 1-5, Pipes A-E).  PipeD is the loop
pipe connecting the two junction nodes; its flow direction emerges from the
hydraulic solution (SIMPLE finds the pressure field that satisfies both
nodal mass balance and the Kirchhoff pressure law around the loop).

The tutorial reports junction pressures and temperatures, per-pipe mass
flows with loop-pipe direction, and the blended mole-fraction composition
at both junctions and the outlet, computed analytically from the converged
flow field.

## [1.6.33] — 2026-08-21

### Tutorials: two compositional solver tutorials; CI fix

**New tutorials** (`tutorials/steady_compositional/`):

- `01_single_pipe/run.py` — methane/ethane binary mixture, 10 km pipe,
  100 → 20 bar.  Prints a PT flash preview (ρ, μ, Cp, k) at inlet and
  outlet conditions, then solves and reports flows, temperatures, and
  gas velocity expansion.
- `02_gas_mixing_junction/run.py` — three-component system (CH₄/C₂H₆/C₃H₈),
  two source nodes at different pressures and compositions, one junction,
  one delivery pipe.  Reports the molar-flow-weighted blend composition
  delivered to the separator.

**CI fix** (`test.yml`): the test workflow now installs `.[dev,compositional]`
so `thermo` is available during CI runs.  The `TestCompositionalFluidProperties`
and `TestSteadyCompositionalSolver` classes also received `@requires_thermo`
skip markers for developers running tests without `thermo` installed.

## [1.6.32] — 2026-08-21

### Feature: compositional pipe network solver (`SteadyCompositionalSolver`)

First fully compositional solver, backed by the `thermo` library (v0.5+) for
equation-of-state PT flashes at each pipe's local pressure and temperature.

**New public API:**

| Symbol | Location |
|--------|----------|
| `CompositionalFluid` | `angelica.properties.compositional_fluid` |
| `CompositionalSolverSettings` | `angelica.solvers.steady_compositional` |
| `SteadyCompositionalSolver` | `angelica.solvers.steady_compositional` |
| `InletCompositionBC` | `angelica.core.case` |

**`CompositionalFluid`** — `FluidModel` implementation that calls
`thermo.Mixture` for PT flashes.  Two-phase VLE regions use the no-slip
homogeneous model (volumetric-fraction-weighted mixture properties).  Results
are cached by `(component_names, P, T, z)` via `functools.lru_cache` to
avoid redundant flash evaluations.

**`SteadyCompositionalSolver`** — Outer loop: propagate compositions →
SIMPLE hydraulics → energy equation → temperature update.  Converges when
`max|Δρ/ρ| < density_rel_tolerance` and `max|ΔT| < temperature_tolerance_k`.

**Multi-inlet composition propagation** — `InletCompositionBC` sets a
mole-fraction Dirichlet condition at each inlet.  In each outer iteration
the solver propagates compositions downstream by BFS-like fixed-point
iteration, mixing streams at junctions by molar-flow-weighted average and
writing the result to `PipeState.zs`.  Handles T-junctions, loops, and
flow reversals automatically.

**Other changes:**

- `PipeState` gains a new field `zs: tuple[float, ...] = ()`.
- `NetworkCase` gains a new field `inlet_composition_bcs: tuple[InletCompositionBC, ...]`.
- New case builder: `angelica.cases.gas_mixing_junction.build_gas_mixing_junction_case`.
- `thermo>=0.5` added as optional dependency (`pip install angelica[compositional]`).
- 24 new tests (all pass); total suite: 199 tests.

## [1.6.31] — 2026-08-21

### Fix: black-oil API now fully exported from top-level package

`SteadyBlackOilSolver`, `BlackOilSolverSettings`, `BlackOilFluid`, `InletFluidBC`,
and `build_water_thermal_fluid` were missing from `angelica/__init__.py`, making
them the only solver type not importable via `from angelica import ...`.

All five are now exported at the top level and listed in `__all__`. The
`__all__` list has also been sorted alphabetically for consistency.

## [1.6.30] — 2026-08-21

### Refactor: remove legacy `examples/` package

Deleted `src/angelica/examples/` (7 thin wrapper modules that re-exported
case builders from `cases/` and called `main()`). The logic has been inlined
directly into the 7 affected tutorial `run.py` files:

- `steady_isothermal_incompressible/01_pipe_only/run.py`
- `steady_isothermal_incompressible/02_fittings_no_elevation/run.py`
- `steady_isothermal_incompressible/03_fittings_with_elevation/run.py`
- `steady_isothermal_incompressible/04_inlet_flow/run.py`
- `steady_isothermal_incompressible/05_outlet_flow/run.py`
- `steady_isothermal_incompressible/06_inlet_and_outlet_flow/run.py`
- `steady_isothermal_incompressible/09_three_reservoir_junction/run.py`

Each tutorial `run.py` now directly imports from `angelica.cases` and
configures its own solver, matching the style of the other tutorials.
All 175 tests pass.

## [1.6.29] — 2026-08-20

### Tests: Hanoi EPANET Darcy-Weisbach benchmark now covered in CI

Added `HanoiBenchmarkTests` to `test_tutorial_suite.py`, covering the last
previously untested tutorial (isothermal incompressible T07).

The test loads the 34-pipe, 32-node Hanoi network from its GUI scene file,
solves it with the settings stored in that file, and asserts:
- All 32 nodal hydraulic heads match the EPANET reference within ±0.05 m.
- Trunk link flows 1–34 match the EPANET reference within ±1 m³/h.

All tutorials across all four solver types now have automated CI coverage.
Total test count: 174 → 175.

## [1.6.28] — 2026-08-20

### Tests: CI coverage for compressible and black-oil tutorial series

Extended `test_tutorial_suite.py` with 12 new tests covering all tutorials
that were previously untested in automated CI.

**`CompressibleTutorialTests` (4 tests)**
- T01 — natural gas pipeline: nodal pressures, pipe flows, and mass balance
- T02 — flow-BC cross-validation: confirms outlet pressures match pressure-BC case
- T03 — looped heat-loss network: convergence, temperature monotone decrease,
  feeder/collector mass balance
- T04 — hill crossing: gravity pressure recovery (P_node4 > P_node2), outflow BC satisfied

**`BlackOilTutorialTests` (6 tests)**
- T01 — single three-phase pipe: flow rate, free-gas appearance below bubble point
- T02 — looped gathering network: loop split, mass balance at T-junctions
- T03 — flow-outlet demand: total throughput equals prescribed 100 kg/s BC
- T04 — two separators: both flow BCs satisfied simultaneously
- T05 — two-reservoir blending: both inlets contribute, mass balance at junction
- T06 — elevation gathering: Well A (lower P, downhill) out-flows Well B (higher P, uphill)

**`NonIsothermalExtendedTests` (2 tests)**
- T05 — crude-oil pipeline thermal: temperature monotone decrease, trunk throughput
- T06 — hilly hot-water network: gravity-assist pressure (P_node3 > P_node1),
  outflow BC satisfied (branch_a = 3.0 kg/s)

Total test count: 162 → 174.

## [1.6.27] — 2026-08-20

### Tutorials: elevation + outflow BC coverage for all solvers

Three new tutorials demonstrate that the elevation source term ρ·g·Δz is
correctly handled across all solver types, including the interaction with
outflow (prescribed-flow) boundary conditions and mixed elevation directions
(some pipes ascending, others descending in the same network).

**Non-isothermal incompressible — Tutorial 06: Hilly hot-water network**
- Branched district-heating network: supply trunk ascends +40 m to a hilltop
  junction, then splits into two downhill branches (−60 m and −30 m).
- Outflow BC at the lower consumer node (ṁ_out = 3.0 kg/s, pressure free).
- Key result: the lower consumer node (−20 m from supply) reaches 694 kPa —
  higher than the 600 kPa inlet — because the 60 m descent from the hilltop
  provides 574 kPa of gravity-assisted driving pressure.

**Compressible — Tutorial 04: Gas pipeline hill crossing**
- Pipeline climbs +500 m to a hilltop junction, then forks into a partial
  descent (−300 m to a pressure outlet) and a full descent past the inlet
  elevation (−550 m to an outflow BC at 8.0 kg/s).
- Elevation correction ρ·g·Δz is re-evaluated each outer iteration as gas
  density varies with pressure.

**Black-oil — Tutorial 06: Gathering network with two wells at different elevations**
- Well A (+150 m, P = 6 MPa) flows downhill to a central manifold — gravity
  provides ~1 MPa of extra driving pressure.
- Well B (−100 m, P = 9 MPa) flows uphill — gravity opposes with ~0.69 MPa.
- Despite its lower wellhead pressure, Well A contributes 58.4 % of total
  production because the elevation advantage more than compensates.

### Improvement: Pipe geometry validation

`Pipe.__post_init__` now raises `ValueError` when
`abs(height_change_m) > length_m`, which is geometrically impossible (a pipe
cannot have a larger elevation change than its own length).  Note that
`length_m` is always the actual pipe axis length (used for friction), not the
horizontal projection — inclined pipes are correctly specified by providing
both `length_m` and `height_change_m` independently.

## [1.6.26] — 2026-08-20

### Performance: sparse pressure solver — O(N³) → O(N·k²) scalability

Replaced `np.linalg.solve` (dense) with `scipy.sparse.linalg.spsolve` (sparse)
for the hydraulic pressure-correction system.

The pressure matrix is a network Laplacian — each node only connects to its
pipe neighbours, so the matrix has at most `2·E` off-diagonal non-zeros (E =
number of pipes). The dense solver allocated O(N²) memory and ran in O(N³)
time; the sparse solver needs only O(N + E) memory and O(N·k²) time where k
is the average node degree (typically 2–4 in pipe networks).

Impact: networks of 10 000+ nodes are now tractable. The EPANET Hanoi
benchmark (34 pipes, 32 nodes) results are numerically unchanged (all 162
tests pass). SciPy was already a required dependency.

## [1.6.25] — 2026-08-20

### Fix: compressible Tutorial 01 crash + reporting test coverage

- `build_natural_gas_pipeline_case()` was missing `thermal_inlets`, causing
  `SteadyCompressibleSolver` to raise `ValueError` ("requires at least one
  ThermalBoundary"). Added `ThermalBoundary(node_id=1, temperature_c=15.0)`,
  consistent with the docstring ("ideal gas, T = 15 °C").
- `test_reporting.py`: renamed the existing test to `test_export_solve_result_workbook_no_balance`
  (tests path with `global_balance=None`, 2-sheet workbook) and added
  `test_export_solve_result_workbook_with_balance_sheet` which provides real
  `GlobalBalance` and `GlobalEnergyBalance` objects and asserts the "Balance"
  sheet is created with the correct row labels.
- Test count: 161 → 162.

## [1.6.24] — 2026-08-20

### Fix: openpyxl bundled with the GUI install

`openpyxl` was an optional extra (`angelica[excel]`) that most users never
installed, causing Excel export to fail with an ImportError. It is now part
of the `[gui]` extras so `pip install 'angelica[gui]'` includes it automatically.
The separate `[excel]` extra is removed.

## [1.6.23] — 2026-08-20

### Feature: global balance included in exported reports

Both CSV and Excel reports now include a **Global Balance** section with the
final mass and energy balance quantities:

- **CSV**: appended as a section at the end of the file.
- **Excel**: written as a dedicated "Balance" sheet.

Quantities reported: mass flow in/out, mass balance error (%); and for thermal
solvers: enthalpy in/out, heat sources, heat wall loss, energy balance error
(kW and %).

## [1.6.22] — 2026-08-20

### Fix: energy balance error now reported as % (consistent with mass balance)

`energy_balance_history` now stores `energy_error_pct` instead of
`abs(energy_error_kw)`, making both balance metrics dimensionless and directly
comparable. The convergence dropdown label is updated to "Energy balance error (%)".

## [1.6.21] — 2026-08-20

### GUI: clean up convergence window and remove balance text panel

- Removed the "Balance" text panel from the sidebar (outside the convergence
  window) — balance information is now exclusively visible as convergence curves.
- Removed the text summary at the bottom of the Convergence Metrics window.
- Removed duplicate mass-balance metrics from the metric dropdown:
  `global_mass_imbalance_kg_per_s` and `global_mass_imbalance_rel` are no
  longer listed (they overlapped with the new per-iteration balance curves).
- Final dropdown now has 6 clean metrics:
  Max |ΔP|, Max nodal mass imbalance, Max |ΔT|, Max |Δρ/ρ|,
  Mass balance error (%), Energy balance error (kW).

## [1.6.20] — 2026-08-20

### Feature: mass and energy balance convergence history in Convergence window

The mass balance error (%) and energy balance error (kW) are now tracked at
every iteration and selectable as metrics in the Convergence Metrics window,
alongside the existing hydraulic metrics.

- **Mass balance error (%)** — `|ṁ_in − ṁ_out| / max(ṁ_in, ṁ_out) × 100`.
  For isothermal runs: tracked at every turbulent inner iteration.
  For non-isothermal / compressible / black-oil runs: tracked at every outer
  (temperature) iteration.
- **Energy balance error (kW)** — `|Ė_in + Q̇_src − Q̇_wall − Ė_out|` in kW.
  Tracked at every outer iteration; zero / not shown for isothermal runs.

Both are exposed as new `mass_balance_history` and `energy_balance_history`
list fields on `SolveResult`.

## [1.6.19] — 2026-08-20

### GUI: energy balance visible in convergence window

The balance panel (mass + energy) is now shown at the bottom of the
Convergence Metrics window, below the plot, in addition to the sidebar.
It is updated automatically when the solve finishes.

## [1.6.18] — 2026-08-20

### Feature: proper global energy balance (first law of thermodynamics)

The global energy balance now tracks the three independent quantities that must
satisfy the steady-state first law for the whole network:

    Ė_in  +  Q̇_src  −  Q̇_wall  =  Ė_out   (residual ≈ 0)

- **Ė_in** — enthalpy rate entering through inlet boundary nodes: Σ ṁ·cp·T
- **Ė_out** — enthalpy rate leaving through outlet boundary nodes: Σ ṁ·cp·T
- **Q̇_src** — heat added by `HeatSource` links (Σ `power_w`)
- **Q̇_wall** — heat lost through pipe walls, computed rigorously:
  - For `n_thermal_segments > 1`: U·π·D·L·(n_internal/n_segs)·(T̄_FV − T_amb),
    matching the FV source terms exactly.
  - For `n_thermal_segments = 1` (NTU analytical bypass): uses the
    log-mean temperature difference (LMTD), which is exact for the bypass path.
- **energy_error_kw** property: residual of the balance; validated at < 1e-8 kW
  for all built-in thermal cases.

A new `GlobalEnergyBalance` dataclass is added to `SolveResult`
(`global_energy_balance`). It is populated by the non-isothermal,
compressible, and black-oil solvers; `None` for the isothermal solver.

The old (incorrect) `heat_loss_kw` field has been removed from `GlobalBalance`,
which now contains only the mass balance.

A "Balance" panel in the GUI sidebar shows the mass balance for all solvers
and the energy balance for thermal solvers, updated live after each run.

## [1.6.17] — 2026-08-20

### Fix: composition convergence checks all four black-oil parameters

- `_propagate_compositions()` previously declared convergence after checking
  only GOR. If two inlets share the same GOR but differ in API gravity, gas
  gravity, or WOR the loop would exit early with wrong downstream compositions.
- Fix: convergence now requires all four parameters (API, gas gravity, GOR,
  WOR) to be stable across all junction nodes.

### Fix: global energy balance uses per-pipe fluid in multi-inlet black-oil networks

- `SteadyBlackOilSolver` was passing `case.fluid_model` (the global reference
  fluid) to `_compute_global_balance`, so the reported heat loss was wrong
  whenever different inlets had different compositions.
- Fix: the call now passes `effective_fluid`, the per-pipe composition-aware
  proxy built in the same outer iteration, so each link's specific heat is
  evaluated at the correct composition.

### Fix: solver runs in background thread — GUI no longer freezes

- The hydraulic solver previously ran on the main thread, blocking the Tkinter
  event loop and making the application unresponsive during long solves.
- Fix: `_run_simulation()` now launches the solver in a `daemon` background
  thread and schedules all UI updates back onto the main thread via
  `root.after(0, ...)`. The convergence plot continues updating live during
  the solve and the window remains interactive.

### CI: test matrix expanded to Python 3.9 – 3.12

- The test workflow previously only ran on Python 3.11. It now runs a parallel
  matrix across 3.9, 3.10, 3.11, and 3.12 with `fail-fast: false`, matching
  the declared `requires-python = ">=3.8"` support range.

## [1.6.16] — 2026-08-20

### Fix: secondary axis curve drawn outside plot area

- The right-axis (green ΔT / Δρ) curve could extend above or below the plot
  boundary when running with tight tolerances (many outer iterations).
- Root cause: the secondary y-axis range was the exact data min/max with no
  margin, so the first (largest) point landed exactly at `y = top` and its
  diamond marker (±4 px) poked outside the canvas.
- Fix: secondary axis now rounds to whole-decade boundaries the same way the
  primary axis does, giving one decade of headroom. `y` is also clamped to
  `[top, bottom]` before drawing, and markers whose centre falls outside the
  plot area are skipped.

## [1.6.15] — 2026-08-19

### GUI: balance residuals selectable in convergence plot

- Replaced the separate balance diagram introduced in v1.6.14 with selectable
  convergence metrics in the existing convergence plot.
- Added per-iteration global mass-balance error metrics:
  `global_mass_imbalance_kg_per_s` and `global_mass_imbalance_rel`.
- Added selectable outer-loop curves for thermal and density convergence:
  `max |ΔT|` and `max |Δρ/ρ|`.

## [1.6.14] — 2026-08-19

### GUI: balance diagram moved to convergence window

- Moved the global mass/energy balance out of the sidebar and into the
  convergence window, where it is shown as a separate diagram next to the
  convergence history.
- The balance view now shows mass entering the network, mass leaving the
  network, net mass error, relative mass error, and thermal energy balance
  when available.
- Updated release metadata so the Windows installer generated from the tag
  reports the current Angelica version.

## [1.6.13] — 2026-08-19

### New: Global mass and energy balance

- After every simulation, a **Global Balance** panel appears in the sidebar
  showing:
  - Total mass flow entering the network (kg/s)
  - Total mass flow leaving the network (kg/s)
  - Relative mass imbalance (%) — should be ~0 for a converged solution
  - Total heat lost to surroundings (kW) — thermal solvers only
    (non-isothermal, compressible, black-oil)
- `SolveResult` gains a `global_balance` field (type `GlobalBalance`) with
  `mass_inlet_kg_per_s`, `mass_outlet_kg_per_s`, `mass_error_pct`, and
  `heat_loss_kw` (None for isothermal).
- Balance is computed in `BaseSolver._compute_global_balance()` from the
  converged network state and is available to all downstream users of
  `SolveResult`.

## [1.6.12] — 2026-08-19

### Fix: title bar now shows the correct version

- `__version__` in `__init__.py` was hardcoded to `"1.3.31"` and never
  updated. It now reads from the installed package metadata via
  `importlib.metadata`, so the title bar always reflects the actual
  installed version.

## [1.6.11] — 2026-08-19

### GUI: Material summary panel simplified for black-oil

- The sidebar material summary now shows only `X°API  γg=Y` in black-oil
  mode. GOR and WOR are per-well properties visible in each source node's
  dialog — they do not belong in the material panel.

## [1.6.10] — 2026-08-19

### GUI: Material dialog simplified to fluid characterization only

- In black-oil mode, Material → Define Material now shows only the two
  fluid characterization fields: **API Gravity** and **Gas Gravity**.
  These are the properties needed to evaluate PVT correlations.
- GOR and WOR are production ratios defined per well — they belong
  exclusively in the source node properties dialog (double-click a
  source node). They no longer appear in the Material dialog at all.

## [1.6.9] — 2026-08-19

### GUI: API gravity and gas gravity moved to Material dialog

- In black-oil mode, **API gravity** and **gas gravity** are now global fluid
  properties set once in Material → Define Material, not per source node.
  These characterize the fluid type and apply to the whole network.
- **GOR** and **WOR** remain per-inlet (source node dialog), since they
  represent individual well production ratios that vary across reservoirs.
- The Material → Define Material dialog for black-oil now shows a "Fluid
  properties" section (API gravity, gas gravity) followed by a
  "Production ratios per inlet" section (GOR and WOR per source node).
- The source node properties dialog now shows only "Production ratios"
  (GOR and WOR), with no fluid characterization fields.
- The material summary panel shows `X°API  γg=Y` on the first line,
  followed by per-inlet GOR/WOR rows.
- Tutorials 01–05 updated: `api_gravity` and `gas_gravity` moved from
  node properties to the `material` dict in each `.gui.json` file.
  Tutorial 05 (two-reservoir blending) now uses a single global fluid
  (32°API, γg=0.65) with different GOR and WOR per reservoir.

## [1.6.8] — 2026-08-19

### GUI: Fluid composition fields restored to source node dialog

- In black-oil mode, double-clicking a source node now shows the four fluid
  composition fields (API gravity, gas gravity, GOR, WOR) directly in the
  node properties dialog, under a "Fluid composition" separator.
- Material → Define Material still provides a summary view of all inlets.
- Both dialogs read and write the same underlying node properties, so they
  stay in sync.

## [1.6.7] — 2026-08-19

### Bug fix

- Black-oil validation error now directs the user to
  "Material → Define Material" instead of "node properties", since fluid
  composition (API, GOR, WOR) was moved to the Material dialog in v1.6.3.

## [1.6.6] — 2026-08-19

### Windows installer fix

- Fixed the Inno Setup script: tutorials were referenced at the old repo-root
  path (`tutorials\*`); updated to `src\angelica\tutorials\*` after the v1.6.2
  move. This caused the Windows installer build to fail for v1.6.2–v1.6.5.
- The build workflow now extracts the version from the git tag and passes it
  to Inno Setup automatically, so `angelica.iss` no longer needs a manual
  version update on each release.

## [1.6.5] — 2026-08-19

### Black-oil tutorials

- Tutorial 05 (Two-Reservoir Blending) now has a GUI file
  (`two_reservoir_blending.gui.json`), making it accessible directly from
  the application. It is the only tutorial that demonstrates per-inlet
  composition with two source nodes carrying different API gravity, GOR,
  and WOR values.
- Added README files for tutorials 03, 04, and 05.
- Added automated test that loads and solves all five black-oil GUI tutorials.

## [1.6.4] — 2026-08-18

### Bug fixes

- `update_node_properties` now merges into existing node properties instead
  of replacing them wholesale. Previously, saving the node BC dialog would
  wipe the fluid composition fields set via the Fluid Definition dialog,
  and vice versa.
- Removed dead `gas_gravity`, `gor`, `wor` widgets from the standard Material
  dialog. They were created visible and immediately hidden by `sync_mode_state`,
  causing a brief flash on open and leaving unreachable code in the dialog.

## [1.6.3] — 2026-08-18

### GUI: Black-oil fluid composition moved to Material dialog

- Composition parameters (API gravity, gas gravity, GOR, WOR) are no longer
  part of the source node properties dialog.
- In black-oil mode, "Material → Define Material" now opens a dedicated
  "Fluid Definition" dialog listing one labeled group per source node, each
  with the four composition fields.
- The material summary panel updates to show each source node's composition
  after saving.

## [1.6.2] — 2026-08-18

### Tutorials now installed with the package

- Moved `tutorials/` inside the Python package (`src/angelica/tutorials/`).
  When Angelica is installed via `pip install`, the tutorials folder is now
  available on all platforms — the GUI "Open" dialog opens there by default.
- Previously the tutorials were only accessible in the Windows installer or
  in a local development clone.

## [1.6.1] — 2026-08-18

### GUI: Composition moved from Material dialog to source node properties

- Removed the "Black-oil" option from the Material dialog — fluid composition
  is no longer a global scene property.
- Each source node in a black-oil network now carries its own four-parameter
  composition (API gravity, gas gravity, GOR, WOR) directly in its Properties
  panel.  This makes it impossible to run a multi-inlet network with a
  single, ambiguous global composition.
- `build_network_case_from_scene()` reads composition from each source node
  and creates the corresponding `InletFluidBC` automatically; the solver path
  is unchanged.
- Tutorials 01–04 updated: the `material` section now contains only a display
  name; composition fields appear on the source node.

## [1.6.0] — 2026-08-17

### New: Per-inlet fluid composition (multi-reservoir black-oil)

- `InletFluidBC` — new dataclass in `angelica.core.case`.  Assigns a
  four-parameter black-oil composition (API gravity, gas gravity, GOR, WOR)
  to a specific inlet node.  Pass one per source node via the new
  `NetworkCase.inlet_fluid_bcs` field.
- `BlackOilComposition` — lightweight frozen dataclass that carries the four
  composition parameters separately from the PVT machinery.  Provides a
  `mix(other, w_self, w_other)` method for mass-weighted blending.
- `compute_pvt()` — standalone module-level function extracted from
  `BlackOilFluid.pvt()`.  Accepts explicit composition parameters so the
  solver can evaluate PVT for any pipe with any composition without needing
  a per-pipe `BlackOilFluid` instance.
- `SteadyBlackOilSolver` now propagates compositions from inlet nodes through
  the network following the flow field, mixes at junctions by mass-weighted
  average, and evaluates density/viscosity per pipe from the local composition.
  Single-fluid networks (no `inlet_fluid_bcs`) behave exactly as before.
- Tutorial 05 — two-reservoir blending: 32°API light crude (Res A, 9 MPa,
  GOR=25) and 22°API heavy crude (Res B, 8 MPa, GOR=10) converge at a mixing
  junction and deliver a blended stream (28.1°API, GOR=19.2) to a separator.
  Confirms that oil, gas, and water conservation holds across the junction.

## [1.5.2] — 2026-08-17

### Tutorials
- Redesign black-oil tutorials 02, 03, and 04 with proper trunk pipes.
  Previously the source node connected directly to the loop T-junction and
  the sink received directly from the T-merge, making the topology look like
  a bare split rather than a realistic gathering network.  All three tutorials
  now include an inlet header (2 km, D=0.22 m) between the source and the
  bifurcation T-junction, and a discharge header (2 km, D=0.22 m) between
  the convergence T-junction and the outlet.  Networks grow from 4 to 6 nodes
  (tutorials 02/03) and from 5 to 7 nodes (tutorial 04).
- Tutorial 04 lateral branch now taps from Node 4 (the lower loop junction),
  which sits between the two trunk pipes, giving the satellite separator a
  physically meaningful upstream pressure (6.10 MPa) above the bubble point.

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
