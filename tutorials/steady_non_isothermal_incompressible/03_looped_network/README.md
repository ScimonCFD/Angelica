# 03 — Looped network with heat loss

## Network

```
source (node 1, 6 bar, 95 °C)
    |
[feed: D=50 mm, L=50 m]
    |
junction A (node 2)
  /                      \
[long path]          [bypass]
2→3: D=20mm, 600m    2→4: D=25mm, 100m
3→4: D=20mm, 400m
  \                      /
junction C (node 4)   ← cold long path + hot bypass mix here
    |
[exit: D=50 mm, L=50 m]
    |
sink (node 5, 1 bar)
```

All pipes: U = 50 W/m²K, T_ambient = 5 °C.

## Why this case requires many outer iterations

The long path (1000 m, D = 20 mm) loses almost all its heat to the
5 °C ambient — junction B arrives near-cold (~21 °C).
The bypass (100 m, D = 25 mm) loses very little heat and stays hot (~91 °C).

This creates a large viscosity contrast between the two paths:
- Long path average ~55 °C → μ ≈ 0.50 mPa·s
- Bypass average ~83 °C → μ ≈ 0.34 mPa·s

The initial temperature guess (all nodes at 95 °C) gives μ ≈ 0.30 mPa·s
everywhere — a completely wrong hydraulic starting point.  As temperatures
update towards their converged values, the flow distribution shifts
significantly, which in turn changes the temperature field.
This mutual coupling keeps the outer loop iterating.

With `temperature_relaxation = 0.5`, the solver takes **14 outer iterations**
with a clean geometric convergence (ratio ≈ 0.5 per step):

```
iter  1: 36.5 K  ████████████████████████████████████████
iter  2: 18.3 K  ████████████████████
iter  3:  9.2 K  ██████████
...
iter 14:  0.005 K
```

## Expected results

| Node | Pressure | Temperature |
|---|---|---|
| Source (1) | 6.000 bar | 95.0 °C |
| Junction A (2) | 5.908 bar | 90.9 °C |
| Junction B — cold end (3) | 3.032 bar | 21.5 °C |
| Junction C — mixing (4) | 1.092 bar | 75.5 °C |
| Sink (5) | 1.000 bar | 72.3 °C |

| Pipe | Flow |
|---|---|
| Feed 1→2 | 1.80 kg/s |
| Long path 2→3 | 0.26 kg/s |
| Long path 3→4 | 0.26 kg/s |
| Bypass 2→4 | 1.54 kg/s |
| Exit 4→5 | 1.80 kg/s |

The bypass carries ~6× more flow than the long path because its
lower viscosity (hotter fluid) and shorter length both reduce resistance.

## Files

| File | Purpose |
|---|---|
| `run.py` | Python script — solves with relax=0.5 and prints convergence bar chart |
| `looped_network_heat_loss.gui.json` | Scene file — open in Angelica GUI |
