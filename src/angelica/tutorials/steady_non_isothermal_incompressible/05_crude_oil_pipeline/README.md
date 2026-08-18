# 05 — Non-isothermal crude oil gathering pipeline

Branched crude oil gathering network demonstrating temperature-dependent fluid properties and viscosity-driven flow behaviour.

**Fluid:** 32°API dead crude oil (Beggs & Robinson viscosity correlation)  
**Network:** one pressure inlet, interior junction, two pressure outlets  
**Inlet:** 800 kPa · 80 °C &nbsp;·&nbsp; **Outlets:** atmospheric (101 325 Pa)  
**Heat loss:** U = 2 W/(m²·K), T_amb = 15 °C

```
inlet (node 1, 800 kPa, 80°C)
    |
[trunk: D=100mm, L=2000m]
    |
junction (node 2)
  /              \
[branch A]      [branch B]
2→3: D=80mm     2→4: D=60mm
     L=1000m         L=600m
     |                |
  outlet A         outlet B
(node 3, atm)   (node 4, atm)
```

Pipes: carbon steel (ε = 0.046 mm), 30 finite-volume thermal segments each.

## What this shows

As crude oil cools from 80 °C toward 15 °C, its viscosity rises by a factor of **44×** (2.9 → 130 cP). A constant-property isothermal solver would severely underestimate friction losses. The non-isothermal solver captures this by updating fluid properties at each cell temperature within the outer temperature iteration.

## Expected output

```
Fluid: 32.0°API dead crude oil
  T =  80.0 °C  →  μ =   2.92 cP  cp = 2103 J/kg/K  k = 0.1205 W/m/K
  T =  15.0 °C  →  μ = 129.65 cP  cp = 1866 J/kg/K  k = 0.1250 W/m/K

Case: Crude Oil Pipeline — thermal (32°API, T_in=80°C)
Converged:              True
Temperature iterations: 4

Node     P (bar)    T (°C)
  1        8.000     80.00
  2        3.716     76.65
  3        1.013     74.36
  4        1.013     74.96

Component         ṁ (kg/s)    Q (m³/h)   T_in (°C)  T_out (°C)
  Pipe:trunk        9.9227     41.3142       80.00       76.65
  Pipe:branch_a      6.1560     25.6310       76.65       74.36
  Pipe:branch_b      3.7667     15.6832       76.65       74.96

Viscosity ratio (cold/hot): 44.4×
```

## Running

```bash
python tutorials/steady_non_isothermal_incompressible/05_crude_oil_pipeline/run.py
```
