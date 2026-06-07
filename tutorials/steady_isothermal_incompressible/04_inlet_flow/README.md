# Tutorial 04 — Inlet Flow Boundary

Same network as tutorial 02 (fittings, no elevation), but node 5 uses a
**mass-flow inlet** instead of a pressure boundary. The mass flow rate equals
the converged value from tutorial 02, so both cases produce the same solution
and can serve as a consistency check.

## Network layout

Identical to tutorial 02.

| Node | Type   | Condition                   |
|------|--------|-----------------------------|
| 1    | Source | 251 300 Pa                  |
| 5    | Source | 0.8477 kg/s (flow inlet)    |
| 6    | Sink   | 101 300 Pa                  |

## What this tutorial demonstrates

- Mass-flow inlet boundary condition
- Mixed boundary types on the same network (pressure + flow)
- Consistency: same converged pressures and flows as tutorial 02

## Run

```bash
python3 tutorials/steady_isothermal_incompressible/04_inlet_flow/run.py
```
