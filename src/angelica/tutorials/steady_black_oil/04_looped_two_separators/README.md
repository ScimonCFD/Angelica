# Tutorial 04 — Black-Oil Loop with Two Separators

Demonstrates a looped black-oil gathering network that splits flow between
**two separator outlets**, each with its own prescribed mass flow rate.

## Network topology

```
N1 (source, P=8 MPa) ──L1(2 km, D=22 cm)──► N2 (T-split)
                                               │
                              ┌────────────────┴───────────────────┐
                              │                                     │
                    L2 (5 km, D=20 cm)                   L4 (5 km, D=18 cm)
                              │                                     │
                              ▼                                     ▼
                           N3 (junction)                       N4 (junction)──L7(2 km, D=14 cm)──► N7 (Separator 2, ṁ=20 kg/s)
                              │                                     │
                    L3 (5 km, D=18 cm)                   L5 (5 km, D=16 cm)
                              │                                     │
                              └────────────────┬───────────────────┘
                                               │
                                           N5 (T-merge)
                                               │
                                    L6 (2 km, D=22 cm)
                                               │
                                               ▼
                                    N6 (Separator 1, ṁ = 80 kg/s)
```

## Boundary conditions

| Node | Type   | Value                |
|------|--------|----------------------|
| N1   | Source | P = 8 MPa, T = 60 °C |
| N6   | Sink   | ṁ = 80 kg/s          |
| N7   | Sink   | ṁ = 20 kg/s          |

## Fluid composition (N1)

| API gravity | Gas gravity | GOR (m³/m³) | WOR (m³/m³) |
|-------------|-------------|-------------|-------------|
| 32 °API     | 0.65        | 25.0        | 0.5         |

## What this tutorial shows

- A single source feeding **two flow-rate sinks** through a looped network.
  The solver finds the inlet pressure and internal flow distribution consistent
  with both outlet flow rates simultaneously.
- The side-branch from N4 to Separator 2 (N7) acts as an additional draw-off
  point mid-loop, making this a more realistic gathering-network topology.

## Running the script

```bash
python run.py
```
