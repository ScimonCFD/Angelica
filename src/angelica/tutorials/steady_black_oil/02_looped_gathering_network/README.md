# Tutorial: Black-Oil Looped Gathering Network

A trunk-fed diamond loop.  A single header feeds a T-split junction, two
parallel loop branches reconnect at a T-merge junction, and a discharge header
leads to the separator.  The solver distributes flow between the branches
according to their hydraulic resistances, and the intermediate pressures
determine the PVT regime (undersaturated vs. two-phase) in each pipe.

## Network

```
                    PipeB  5 km, D=0.20 m          PipeC  5 km, D=0.18 m
              ┌────────────────────── Node 3 ──────────────────────┐
              │                                                     │
  Node 1 ─PipeA─ Node 2                                         Node 5 ─PipeF─ Node 6
  8 MPa   2 km   T-split                                        T-merge  2 km  2 MPa
  60 °C   0.22 m  │                                                     │  0.22 m
              └────────────────────── Node 4 ──────────────────────┘
                    PipeD  8 km, D=0.15 m          PipeE  3 km, D=0.20 m
```

| Node | Role | Boundary |
|---|---|---|
| 1 | Inlet | P = 8 MPa, T = 60 °C |
| 2 | T-split | Free |
| 3 | Upper loop junction | Free |
| 4 | Lower loop junction | Free |
| 5 | T-merge | Free |
| 6 | Separator outlet | P = 2 MPa |

Upper branch (2 → 3 → 5): 10 km, wider pipes
Lower branch (2 → 4 → 5): 11 km, narrow 8 km section (PipeD)

## Fluid

32 °API crude, γ_g = 0.65, GOR = 25 m³/m³, WOR = 0.5 m³/m³.  
Bubble point at 60 °C: **5.63 MPa**.

## Key results

| Node | P (MPa) | Regime |
|---|---|---|
| 1 | 8.000 | Undersaturated |
| 2 | 7.176 | Undersaturated |
| 3 | 5.793 | Undersaturated (just above Pb) |
| 4 | 3.567 | Two-phase |
| 5 | 3.174 | Two-phase |
| 6 | 2.000 | Two-phase |

| Branch | Mass flow | Share |
|---|---|---|
| Upper (B + C) | 75.7 kg/s | 62.9 % |
| Lower (D + E) | 44.7 kg/s | 37.1 % |

**PVT effect on volumetric flow:** Pipe B and Pipe C carry the same mass but
different mixture volumes (311 vs 343 m³/h).  Pipe B stays above the bubble
point (Node 3 is at 5.79 MPa), while along Pipe C the pressure falls below Pb,
liberating free gas and expanding the mixture.  In the lower branch the narrow
PipeD already takes the pressure below Pb (Node 4 at 3.57 MPa), so both lower
pipes carry a two-phase mixture (191 and 225 m³/h).  The discharge header
carries 687 m³/h of mixture for the 495 m³/h entering the trunk.

Solver converges in **3 PVT iterations**.

## Run

```bash
python3 tutorials/steady_black_oil/02_looped_gathering_network/run.py
```
