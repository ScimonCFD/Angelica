# Tutorial 03 — Black-Oil Loop with Flow Outlet BC

Demonstrates a looped black-oil gathering network where the outlet boundary
condition is a **prescribed mass flow rate** instead of a fixed pressure.

## Network topology

```
N1 (source, P=8 MPa) ──L1(2 km, D=22 cm)──► N2 (T-split)
                                               │
                              ┌────────────────┴───────────────────┐
                              │                                     │
                    L2 (5 km, D=20 cm)                   L4 (8 km, D=15 cm)
                              │                                     │
                              ▼                                     ▼
                           N3 (junction)                       N4 (junction)
                              │                                     │
                    L3 (5 km, D=18 cm)                   L5 (3 km, D=20 cm)
                              │                                     │
                              └────────────────┬───────────────────┘
                                               │
                                           N5 (T-merge)
                                               │
                                    L6 (2 km, D=22 cm)
                                               │
                                               ▼
                                    N6 (sink, ṁ = 100 kg/s)
```

## Boundary conditions

| Node | Type   | Value             |
|------|--------|-------------------|
| N1   | Source | P = 8 MPa, T = 60 °C |
| N6   | Sink   | ṁ = 100 kg/s      |

## Fluid composition (N1)

| API gravity | Gas gravity | GOR (m³/m³) | WOR (m³/m³) |
|-------------|-------------|-------------|-------------|
| 32 °API     | 0.65        | 25.0        | 0.5         |

## What this tutorial shows

- How to use a **flow-rate outlet** (instead of a pressure outlet) in a
  black-oil looped network.  The solver finds the inlet pressure consistent
  with delivering the specified total flow.
- The loop naturally distributes the flow between the upper path (N2→N3→N5)
  and the lower path (N2→N4→N5) according to the hydraulic resistance of each
  branch.
- Trunk pipes before and after the loop make the topology explicit.

## Running the script

```bash
python run.py
```
