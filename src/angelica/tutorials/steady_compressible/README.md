# steady_compressible

Tutorial cases for the steady, single-phase compressible solver.

The solver uses an outer density loop that wraps the SIMPLE hydraulic solver:
after each full hydraulic solve, fluid density is re-evaluated from the equation
of state at the updated pressure field.  The loop repeats until the maximum
relative density change across all links falls below the convergence tolerance.

**Scope:** single-phase gases described by an equation of state (ideal gas
or Peng-Robinson). Temperature is solved with the flow: wall heat exchange,
compressor shaft work, Joule-Thomson cooling across valves and — since
v1.18.0 — along pipes, and the potential energy of elevation changes. It is
not valid for two-phase flow or near the critical point (use the
compositional solver).

Available cases:

- `01_natural_gas_pipeline` — branched methane network, 800 → 500 kPa; density gradient along the lines.
- `02_flow_bc_cross_validation` — the same network driven by a flow boundary instead of a pressure boundary gives the same solution.
- `03_looped_gas_pipeline_heat_loss` — looped gas network with wall heat loss.
- `04_hill_crossing` — gas over a 500 m hill: the elevation term with a pressure-dependent density.
- `05_valve_jt_effect` — Joule-Thomson cooling across a throttling valve and in the pipes before and after it.
- `06_compressor_heating` — compressor station: shaft-work heating of a real gas; Part B
  compares the head-curve machine with a `Compressor` (isentropic and polytropic
  efficiency, power and discharge temperature from the gas's isentrope).

## Running a tutorial

```bash
python tutorials/steady_compressible/01_natural_gas_pipeline/run.py
```

Or open the `.gui.json` file in Angelica via **File → Open**.
