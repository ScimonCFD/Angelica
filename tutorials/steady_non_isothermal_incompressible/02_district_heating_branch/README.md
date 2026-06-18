# 02 — District heating branch

## Network

```
source (node 1, 4 bar, 85 °C)
    |
[main pipe: D=50 mm, L=200 m, U=5 W/m²K]
    |
junction (node 2)
  /              \
[branch A]      [branch B]
D=25 mm, L=300 m   D=25 mm, L=100 m
U=5 W/m²K          U=5 W/m²K
  |                   |
sink A (node 3)   sink B (node 4)
  1.5 bar             1.5 bar
```

## Parameters

| Property | Value |
|---|---|
| Fluid | Water (ρ=998 kg/m³, μ=0.001 Pa·s, cₚ=4182 J/kg·K, k=0.6 W/m·K) |
| Ambient temperature | 10 °C |
| Supply pressure | 4 bar |
| Delivery pressure | 1.5 bar (both sinks) |

## What to expect

| Location | Pressure | Temperature |
|---|---|---|
| Supply (node 1) | 4.00 bar | 85.0 °C |
| Junction (node 2) | ~3.69 bar | ~83.3 °C |
| Load A — long branch (node 3) | 1.5 bar | ~79.9 °C |
| Load B — short branch (node 4) | 1.5 bar | ~82.7 °C |

The longer branch (300 m) loses more heat than the shorter one (100 m),
so load A receives water at a lower temperature than load B even though
both start from the same junction.

## Files

| File | Purpose |
|---|---|
| `run.py` | Python script — solves and prints a summary table |
| `district_heating_branch.gui.json` | Scene file — open in Angelica GUI |
