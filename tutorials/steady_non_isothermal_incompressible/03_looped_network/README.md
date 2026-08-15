# 03 — Looped network with heat loss

## Network

```
source (node 1, 5 bar, 95 °C)
    |
[feed: D=50 mm, L=100 m]
    |
junction A (node 2)
  /                      \
[upper branch]        [direct bypass]
2→3: D=40mm, 400m     2→4: D=35mm, 400m
3→4: D=40mm, 300m
  \                      /
junction C (node 4)   ← upper branch + bypass mix here
    |
[exit: D=50 mm, L=100 m]
    |
sink (node 5, 1 bar)
```

All pipes: U = 30 W/m²K, T_amb = 10 °C.

## Parameters

| Property | Value |
|---|---|
| Fluid | Water (temperature-dependent ρ, μ) |
| Inlet pressure | 5 bar |
| Outlet pressure | 1 bar |
| Inlet temperature | 95 °C |
| Ambient temperature | 10 °C |
| Heat transfer coefficient U | 30 W/m²K |

## What to expect

Hot water enters at 95 °C and splits at junction A into two paths:
- **Upper branch** (700 m total, D=40 mm): longer, loses more heat
- **Direct bypass** (400 m, D=35 mm): shorter, loses less heat

Both paths rejoin at junction C, where their temperatures mix.
All temperatures decrease monotonically along each flow path.

The two paths have slightly different average temperatures after splitting,
so their fluid viscosities differ.  Since viscosity affects friction and
therefore flow distribution, the hydraulic and thermal solutions are coupled
— the outer temperature loop must iterate several times to converge.

## Expected results

| Node | Pressure | Temperature |
|---|---|---|
| Source (1) | 5.000 bar | 95.0 °C |
| Junction A (2) | 4.507 bar | 91.9 °C |
| Junction B (3) | 2.789 bar | 74.9 °C |
| Junction C (4) | 1.495 bar | 70.0 °C |
| Sink (5) | 1.000 bar | 67.8 °C |

| Pipe | Flow |
|---|---|
| Feed 1→2 | 3.00 kg/s |
| Upper 2→3 | 1.55 kg/s |
| Upper 3→4 | 1.55 kg/s |
| Bypass 2→4 | 1.45 kg/s |
| Exit 4→5 | 3.00 kg/s |

## Files

| File | Purpose |
|---|---|
| `run.py` | Python script — solves and prints results plus convergence history |
| `looped_network_heat_loss.gui.json` | Scene file — open in Angelica GUI |
