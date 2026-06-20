# 01 — Single pipe with heat loss

## Network

```
source (node 1)  ──[pipe, 500 m, D=25 mm]──>  sink (node 2)
  3 bar, 80 °C                                    1 bar
```

## Parameters

| Property | Value |
|---|---|
| Fluid | Water (ρ=998 kg/m³, μ=0.001 Pa·s, cₚ=4182 J/kg·K, k=0.6 W/m·K) |
| Pipe diameter | 25 mm |
| Pipe length | 500 m |
| Wall heat transfer coeff. U | 50 W/m²K |
| Ambient temperature | 20 °C |
| Thermal segments | 50 |

## What to expect

Hot water enters at 80 °C and loses heat through the pipe wall.
At the converged mass flow (~0.41 kg/s), the outlet temperature is roughly
**39–41 °C** depending on mesh resolution.

The analytical solution for plug-flow with a uniform heat sink is:

```
T_out = T_amb + (T_in − T_amb) × exp(−U·π·D·L / (ṁ·cₚ))
```

This gives ~38.8 °C at ṁ=0.41 kg/s; the FV result converges towards this as
the number of thermal segments is increased.

## Files

| File | Purpose |
|---|---|
| `run.py` | Python script — solves and prints results plus analytical comparison |
| `hot_water_pipe_heat_loss.gui.json` | Scene file — open in Angelica GUI |
