# Tutorial 04 — Inline electric heater

Cold water at 20 °C flows through a 50 mm pipe driven by a 2 bar pressure
difference. An inline electric heater adds 50 kW to the fluid. The heater
is hydraulically transparent (ΔP = 0).

**Network layout**

```
Source (3 bar, 20 °C)
  ── Pipe (D=50 mm, L=10 m, insulated) ──
  ── HeatSource (Q = 50 kW, ΔP = 0) ──
  ── Pipe (D=50 mm, L=10 m, insulated) ──
Sink (1 bar)
```

**Expected result**

The temperature rise across the heater satisfies the steady-state energy balance:

```
ΔT = Q / (ṁ · cp)
```

At the flow rate produced by the 2 bar driving pressure (~13.9 kg/s), ΔT ≈ 0.86 K.
The solver matches this to within 0.5 K.

## Running

```bash
python tutorials/steady_non_isothermal_incompressible/04_inline_heater/run.py
```

Or open `inline_heater.gui.json` in Angelica via **File → Open**.
