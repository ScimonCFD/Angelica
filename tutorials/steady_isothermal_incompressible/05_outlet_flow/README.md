# Tutorial 05 — Outlet Flow Boundary

Same network as tutorial 03 (fittings with elevation), but node 6 uses a
**mass-flow outlet** instead of a pressure boundary. The mass flow rate equals
the converged value from tutorial 03, so both cases produce the same solution
and serve as a consistency check.

## Network layout

Identical to tutorial 03 (fittings + elevation changes).

| Node | Type   | Condition                    |
|------|--------|------------------------------|
| 1    | Source | 251 300 Pa                   |
| 5    | Source | 201 300 Pa                   |
| 6    | Sink   | 1.7717 kg/s (flow outlet)    |

## What this tutorial demonstrates

- Mass-flow outlet boundary condition
- Combined use of two pressure inlets and one flow outlet
- Consistency: same converged pressures and flows as tutorial 03

## Run

```bash
python3 tutorials/steady_isothermal_incompressible/05_outlet_flow/run.py
```
