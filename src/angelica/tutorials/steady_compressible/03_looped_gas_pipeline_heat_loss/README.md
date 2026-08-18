# 03 — Looped gas pipeline with heat loss

## Network

```
source (node 1, 700 kPa, 40 °C)
    |
[feeder: D=200 mm, 600 m]
    |
junction A (node 2)
  /                    \
[upper]              [lower]
2→3: D=150mm, 500m   2→4: D=120mm, 400m
       |                     |
      L4: D=100mm, 350m      |
       |                     |
3→5: D=120mm, 400m   4→5: D=150mm, 500m
  \                    /
junction B (node 5)
    |
[collector: D=200 mm, 600 m]
    |
sink (node 6, 500 kPa)
```

Gas: methane (CH₄). All pipes: roughness = 0.046 mm, U = 5 W/m²K, T_amb = 10 °C.

## Parameters

| Property | Value |
|---|---|
| Gas | Methane (CH₄, M = 16.04 g/mol) |
| Inlet pressure | 700 kPa |
| Outlet pressure | 500 kPa |
| Inlet temperature | 40 °C |
| Ambient temperature | 10 °C |
| Heat transfer coefficient U | 5 W/m²K |
| Pipe roughness | 0.046 mm |

## What to expect

Methane enters at 700 kPa and 40 °C through a large feeder pipe and splits at
junction A into two asymmetric paths:

- **Upper path** (2→3→5): 500 m + 400 m, smaller-to-larger diameters
- **Lower path** (2→4→5): 400 m + 500 m, larger-to-smaller diameters

The diagonal cross-pipe (L4: 3→4) allows flow to redistribute between the
two branches, equalising pressure at junctions 3 and 4.  Because the upper
and lower paths are mirror-asymmetric (different length/diameter combinations),
the cross-pipe carries real flow — any perfectly symmetric design would give
zero cross-flow, making the loop redundant.

Gas cools as it flows, with most of the temperature drop occurring in the
long feeder and collector pipes.

## Expected results

| Link | Route | Flow (kg/s) |
|---|---|---|
| L1 | feeder 1→2 | 2.639 |
| L2 | upper 2→3 | 1.518 |
| L3 | lower 2→4 | 1.122 |
| L4 | cross 3→4 | 0.395 |
| L5 | upper 3→5 | 1.123 |
| L6 | lower 4→5 | 1.517 |
| L7 | collector 5→6 | 2.639 |

| Node | Pressure | Temperature |
|---|---|---|
| Source (1) | 700.0 kPa | 40.0 °C |
| Junction A (2) | 663.4 kPa | 32.6 °C |
| Junction 3 | 616.4 kPa | 26.6 °C |
| Junction 4 | 596.9 kPa | 25.2 °C |
| Junction B (5) | 545.3 kPa | 21.8 °C |
| Sink (6) | 500.0 kPa | 18.9 °C |

## Running

```bash
python tutorials/steady_compressible/03_looped_gas_pipeline_heat_loss/run.py
```

Or open `looped_gas_pipeline_heat_loss.gui.json` in Angelica via **File → Open**.

## Files

| File | Purpose |
|---|---|
| `run.py` | Scripted API example |
| `looped_gas_pipeline_heat_loss.gui.json` | Scene file — open in Angelica GUI |
