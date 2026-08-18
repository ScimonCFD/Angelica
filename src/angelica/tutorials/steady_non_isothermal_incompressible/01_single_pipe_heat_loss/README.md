# 01 — Single pipe with heat loss

## Network

```
source (node 1)  ──[pipe, 500 m, D=25 mm]──>  sink (node 2)
  2 bar, 80 °C                                    1 bar
```

## Parameters

| Property | Value |
|---|---|
| Fluid | Water (ρ=998 kg/m³, μ=0.001 Pa·s, cₚ=4182 J/kg·K, k=0.6 W/m·K) |
| Pipe diameter | 25 mm |
| Pipe length | 500 m |
| Wall heat transfer coeff. U | 50 W/m²K |
| Ambient temperature | 20 °C |

## What to expect

Hot water enters at 80 °C and loses heat through the pipe wall.
The analytical solution for plug-flow with a uniform heat sink is:

```
T_out = T_amb + (T_in − T_amb) × exp(−U·π·D·L / (ṁ·cₚ))
```

At the converged mass flow (~0.41 kg/s) this gives **~38.8 °C**.
The solver uses the NTU formula directly for each pipe, so the result
is exact regardless of the number of pipe segments.

## Files

| File | Purpose |
|---|---|
| `run.py` | Python script — solves and prints results plus analytical comparison |
| `hot_water_pipe_heat_loss.gui.json` | Scene file — open in Angelica GUI |
