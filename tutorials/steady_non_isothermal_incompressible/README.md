# steady_non_isothermal_incompressible

Tutorial cases for the steady, incompressible, non-isothermal solver.

The solver uses an outer temperature loop that wraps the SIMPLE hydraulic
solver: hydraulics are solved with the current density and viscosity, then
the energy equation is solved on the resulting flow field, and temperatures
are updated until convergence.

Available cases:

- `01_single_pipe_heat_loss` — single pipe analytical benchmark (NTU method)
- `02_district_heating_branch` — branched network with two thermal loads
- `03_looped_network` — looped network with ambient heat loss and temperature convergence
- `04_inline_heater` — inline electric heater, energy balance verification

## Running a tutorial

```bash
python tutorials/steady_non_isothermal_incompressible/01_single_pipe_heat_loss/run.py
```

Or open the `.gui.json` file directly in Angelica via **File → Open**.
