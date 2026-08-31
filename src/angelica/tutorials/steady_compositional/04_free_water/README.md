# Tutorial 04 — Free Water in a Wet-Gas Pipeline

Demonstrates Angelica's **immiscible-water model**.  When `"water"` (or
`"h2o"`, `"7732-18-5"`) is included in the component list the solver detects
it automatically and applies a two-step approach:

1. VL flash is performed on the **dry** (water-free) normalised composition.
2. The water phase is split using the **Wagner / IAPWS-IF97 saturation
   pressure**: water condenses when its partial pressure in the vapour phase
   would exceed `Psat(T)`.

No API change is required — just add water to the component and mole-fraction
vectors.

## Geometry

```
Node 1 ──── Pipe A (50 km, D = 0.20 m) ──── Node 2 ──── Pipe B (50 km) ──── Node 3
P = 80 bar, T = 50 °C                                                         P = 20 bar
```

## Fluid

| Component | mol% |
|-----------|-----:|
| methane   | 75.0 |
| ethane    |  8.0 |
| propane   |  4.0 |
| n-butane  |  1.0 |
| **water** | **12.0** |

## What the tutorial shows

| Section | Description |
|---------|-------------|
| 1 | PVT comparison: wet gas vs dry gas flash at inlet conditions |
| 2 | Free-water fraction scan across T and P |
| 3 | Full network solve with temperature cooling |
| 4 | `free_water_fraction` per pipe and estimated liquid-water mass flow |

## Key result

At 80 bar the water saturation pressure is only ~0.12 bar, so almost all of
the 12 mol% water feed condenses as free liquid water.  As the gas cools
along the pipeline, slightly more water condenses (Psat decreases with T).

## Running

```bash
python run.py
```
