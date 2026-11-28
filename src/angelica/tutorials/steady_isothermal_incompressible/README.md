# steady_isothermal_incompressible

This folder contains the tutorial suite for the steady, incompressible,
isothermal, single-component solver.  What each case is checked against —
an analytical solution (10 Poiseuille, 09 three reservoirs), a published
benchmark (08 EPANET pump), or only Angelica's own stored output (the
rest, including the Hanoi network) — is listed in the validation manual,
`docs/validation.md`.

The available cases are:

- `01_pipe_only`
- `02_fittings_no_elevation`
- `03_fittings_with_elevation`
- `04_inlet_flow`
- `05_outlet_flow`
- `06_inlet_and_outlet_flow`
- `07_hanoi_epanet_darcy_benchmark`
- `08_epanet_pump_benchmark`
- `09_three_reservoir_junction`
- `10_laminar_poiseuille_benchmark`
- `11_crude_oil_pipeline`
