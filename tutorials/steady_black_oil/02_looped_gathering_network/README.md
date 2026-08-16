# Tutorial: Black-Oil Looped Gathering Network

A looped pipeline network with two parallel flow paths between a high-pressure
inlet and a low-pressure separator.  The solver distributes flow between the
paths according to their hydraulic resistances, and the intermediate pressures
determine the PVT regime (undersaturated vs. two-phase) in each branch.

## Network

```
                PipeA (5 km, D=0.20 m)
  [Node 1] ─────────────────────────────> [Node 2] ──┐
      │                                              │
      │  PipeC (8 km, D=0.15 m)                  PipeB
      └────────────────────────> [Node 3]         (5 km, D=0.18 m)
                                    │                │
                                 PipeD               │
                              (3 km, D=0.20 m)       │
                                    └────────────> [Node 4]
                                                  Separator
```

| Node | Role | Boundary |
|---|---|---|
| 1 | Inlet | P = 8 MPa, T = 60 °C |
| 2 | Upper junction | Free |
| 3 | Lower junction | Free |
| 4 | Separator outlet | P = 2 MPa |

Upper path (1 → 2 → 4): 10 km total, wider diameters  
Lower path (1 → 3 → 4): 11 km total, narrower middle section

## Fluid

32 °API crude, γ_g = 0.65, GOR = 25 m³/m³, WOR = 0.5 m³/m³.  
Bubble point at 60 °C: **5.63 MPa**.

## Key results

| Node | P (MPa) | Regime |
|---|---|---|
| 1 | 8.00 | Undersaturated |
| 2 | 5.99 | Undersaturated (just above Pb) |
| 3 | 2.68 | Two-phase |
| 4 | 2.00 | Two-phase |

| Path | Mass flow | Share |
|---|---|---|
| Upper (A + B) | 92.4 kg/s | 62.7 % |
| Lower (C + D) | 54.9 kg/s | 37.3 % |

**PVT effect on volumetric flow:** Pipe A and Pipe B carry the same mass but
different mixture volumes (380 vs 439 m³/h).  The increase happens because
the pressure in Pipe B falls below the bubble point along the way, liberating
free gas and expanding the mixture.  This is not observed in the lower path
because Pipe C is the narrower bottleneck — the pressure at Node 3 is already
below Pb, so the lower path is two-phase from the start.

Solver converges in **3 PVT iterations**.

## Run

```bash
python3 tutorials/steady_black_oil/02_looped_gathering_network/run.py
```
