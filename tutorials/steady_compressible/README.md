# steady_compressible

Tutorial cases for the steady, single-phase compressible solver.

The solver uses an outer density loop that wraps the SIMPLE hydraulic solver:
after each full hydraulic solve, fluid density is re-evaluated from the equation
of state at the updated pressure field.  The loop repeats until the maximum
relative density change across all links falls below the convergence tolerance.

**Scope:** This solver is designed for **single-phase gases** modelled by an
equation of state — currently ideal gas (ρ = PM/RT).  It is not valid for
two-phase flow or near the critical point.  Temperature is treated as uniform
(isothermal); a non-isothermal compressible solver is planned.

Available cases:

- `01_natural_gas_pipeline` — branched methane network, 800 → 500 kPa, demonstrates density gradient

## Running a tutorial

```bash
python tutorials/steady_compressible/01_natural_gas_pipeline/run.py
```

Or open the `.gui.json` file in Angelica via **File → Open**.
