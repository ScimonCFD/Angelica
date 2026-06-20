# steady_non_isothermal_incompressible

Tutorial cases for the steady, incompressible, non-isothermal solver.

The solver uses an outer temperature loop that wraps the SIMPLE hydraulic
solver: hydraulics are solved with the current density and viscosity, then
the energy equation is solved on the resulting flow field, and temperatures
are updated until convergence.

**Scope:** This solver is designed for **incompressible liquids** (water,
thermal oils, crude oil) where fluid properties depend on temperature but not
on pressure.  It is not valid for compressible gases, where density is a
function of both pressure and temperature.  Transient effects are not
modelled.

Available cases:

- `01_single_pipe_heat_loss` — single pipe analytical benchmark (NTU method)
- `02_district_heating_branch` — branched network with two thermal loads
- `03_looped_network` — looped network with ambient heat loss and temperature convergence
- `04_inline_heater` — inline electric heater, energy balance verification

Additional validation cases currently covered by the automated test suite:

- `thermal_mixing_junction` — exact adiabatic junction mixing benchmark
- `inline_heater_fixed_flow` — exact `ΔT = Q / (ṁ cₚ)` benchmark with prescribed mass flow
- `symmetric_adiabatic_loop` — looped network with exact 50/50 split and uniform temperature
- `symmetric_heat_loss_loop` — looped network with symmetric NTU heat loss and exact mixed outlet temperature

## Running a tutorial

```bash
python tutorials/steady_non_isothermal_incompressible/01_single_pipe_heat_loss/run.py
```

Or open the `.gui.json` file directly in Angelica via **File → Open**.
