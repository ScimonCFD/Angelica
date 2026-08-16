# Tutorial: Black-Oil Three-Phase Pipeline

A single crude-oil pipe transporting a live-oil mixture (dissolved gas + produced water)
from a high-pressure inlet to a low-pressure outlet.  The tutorial illustrates
two classic PVT regimes in the same pipe: **undersaturated** flow near the inlet
(pressure above the bubble point, all gas dissolved) and **two-phase** flow near
the outlet (pressure below the bubble point, free gas appears).

## Network

```
[Node 1] ──── Pipe · 10 km · D = 0.2 m ──── [Node 2]
 P = 8 MPa                                    P = 2 MPa
 T = 60 °C
```

## Fluid

| Property | Value |
|---|---|
| API gravity | 32 °API |
| Gas gravity (air = 1) | 0.65 |
| GOR at standard conditions | 25 m³/m³ |
| WOR at standard conditions | 0.5 m³/m³ |

## PVT behaviour

Bubble point at 60 °C: **≈ 5.6 MPa** (Standing 1947).

| Location | P | Regime | Gas holdup |
|---|---|---|---|
| Inlet | 8 MPa | Undersaturated (P > Pb) | 0 % |
| Outlet | 2 MPa | Two-phase (P < Pb) | ≈ 40 % |

Because the inlet is above the bubble point, the gas holdup there is exactly zero
and the mixture behaves as a two-component (oil + water) liquid.  As pressure
drops along the pipe, gas liberates from solution and the mixture density falls
from ≈ 875 kg/m³ at the inlet to ≈ 540 kg/m³ at the outlet.

## Correlations

| Quantity | Correlation |
|---|---|
| Bubble point, Rs, Bo | Standing (1947) |
| Gas z-factor | Hall-Yarborough (1974) + Sutton (1985) pseudo-crits |
| Gas viscosity | Lee, Gonzalez & Eakin (1966) |
| Live-oil viscosity | Beggs & Robinson (1975) |
| Water FVF | McCain (1990), simplified |

## Run

```bash
python3 tutorials/steady_black_oil/01_three_phase_pipeline/run.py
```

## Expected output (abridged)

```
Bubble point at 60 °C : 5.63 MPa
  → Inlet (8 MPa) is UNDERSATURATED
  → Outlet (2 MPa) is TWO-PHASE

                                 Inlet (8 MPa)  Outlet (2 MPa)
Rs (m³/m³)                               25.00            6.75
Bo (m³_res/m³_sc)                       1.0901          1.0497
Holdup gas                              0.0000          0.3989
ρ_mixture (kg/m³)                        875.2           539.2

Converged: True  |  PVT iterations: 2
Mass flow: 109.79 kg/s  |  T outlet: 56.08 °C

Surface rates:
  Oil:   282.92 m³/h
  Gas:   7072.9 m³/h  (25.0 m³/m³ ← GOR input: 25.0)
  Water: 141.46 m³/h  (0.50 m³/m³ ← WOR input: 0.5)
```

The solver converges in 2 PVT outer iterations.  Surface GOR and WOR recover
exactly from the mass-flow split because the surface composition is fixed by
the input parameters.
