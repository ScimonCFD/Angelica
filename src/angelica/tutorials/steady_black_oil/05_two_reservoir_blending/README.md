# Tutorial 05 — Black-Oil Two-Reservoir Blending

Demonstrates the **per-inlet fluid composition** capability: two reservoirs
with different crude properties feed a common trunk, and the solver propagates
and mixes their compositions at the junction.

## Network topology

```
N1 (Reservoir A, P=9 MPa, T=70 °C)
  │  32°API  GOR=25  WOR=0.3
  │
  └──L1 (5 km, D=18 cm)──► N3 (Mixing point) ──L3 (3 km, D=22 cm)──► N4 (Separator, P=2 MPa)
                             ▲
  ┌──L2 (5 km, D=16 cm)─────┘
  │
N2 (Reservoir B, P=8 MPa, T=60 °C)
  22°API  GOR=10  WOR=1.5
```

## Boundary conditions

| Node | Type   | Value                 |
|------|--------|-----------------------|
| N1   | Source | P = 9 MPa, T = 70 °C  |
| N2   | Source | P = 8 MPa, T = 60 °C  |
| N4   | Sink   | P = 2 MPa             |

## Fluid compositions

| Reservoir | API gravity | Gas gravity | GOR (m³/m³) | WOR (m³/m³) |
|-----------|-------------|-------------|-------------|-------------|
| A (N1)    | 32 °API     | 0.65        | 25.0        | 0.3         |
| B (N2)    | 22 °API     | 0.70        | 10.0        | 1.5         |

## What this tutorial shows

- **Per-inlet composition**: each source node carries its own fluid
  definition.  The solver assigns PVT properties per pipe based on which
  reservoir's fluid occupies it.
- **Composition mixing**: at the junction (N3) the solver computes a
  mass-weighted average of the incoming compositions. The blended fluid
  (API, GOR, WOR) flows through the trunk to the separator.
- The pressure difference between the two reservoirs drives unequal flow
  rates, so the blended composition at N3 reflects the actual mass split,
  not a simple 50/50 average.

## Expected results (approximate)

The heavier reservoir (B, 22°API) delivers more water but less gas per unit
mass than reservoir A. The blended stream at the separator has intermediate
API gravity and GOR, weighted by the respective flow rates.

## Running the script

```bash
python run.py
```
