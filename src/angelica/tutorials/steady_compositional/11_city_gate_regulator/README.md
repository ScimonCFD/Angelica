# Tutorial 11 — A pressure-reducing station (city gate)

A dry natural gas arrives through a 20 km transmission line and is let down
to a 20 bar distribution network by a pressure-reducing regulator.

```
Supply 70 bar, 20 °C ─ transmission 20 km, 8" ─ [heater] ─ Regulator (20 bar) ─ distribution 2 km, 6" ─ Consumers 5 bar
```

Gas (mole fractions): C1 0.90, C2 0.05, C3 0.02, nC4 0.01, N2 0.01,
CO2 0.01.  Ground at 10 °C, U = 2 W/(m²·K).

## How it is modelled

- `PressureRegulator(..., set_pressure_pa=20e5)` (mode "downstream") holds
  its outlet at 20 bar while it can, and passes what the distribution side
  takes.  When the supply cannot keep its inlet high enough it is wide open
  (a fitting, K = 5); when the outlet side is at or above 20 bar it is
  closed.  See black-oil tutorial 15 for the three states.
- The gas crosses it at constant enthalpy (Peng-Robinson): it cools.
- The heater is a fixed-duty `HeatSource` before the regulator.

## Results

**Part A — let-down without preheating.**
- The regulator controls: 63.79 → 20.00 bar, 7.589 kg/s.
- The gas cools from 11.44 to −15.79 °C across it (−0.62 K/bar on average,
  the Peng-Robinson isenthalpic expansion of this gas).
- It reaches the consumers at −23.6 °C: the distribution line expands it
  further.
- Energy balance residual 0.003 kW.

**Part B — preheating.**

| Heater (kW) | T in (°C) | T out (°C) | Flow (kg/s) |
|---|---|---|---|
| 0 | 11.44 | −15.79 | 7.589 |
| 200 | 21.59 | −3.68 | 7.386 |
| 400 | 32.50 | 9.14 | 7.191 |
| 600 | 44.18 | 22.65 | 7.004 |

About 260 kW keeps the regulator outlet at 0 °C.  The warmer gas is less
dense, so the distribution line passes less mass at the same pressures.

**Part C — supply pressure.**

| Supply (bar) | State | P in (bar) | P out (bar) | Flow (kg/s) |
|---|---|---|---|---|
| 25 | open | 14.53 | 14.48 | 4.995 |
| 30 | open | 17.24 | 17.18 | 6.077 |
| 35 | open | 19.98 | 19.90 | 7.169 |
| 40 | controlling | 27.69 | 20.00 | 7.263 |
| 50 | controlling | 40.85 | 20.00 | 7.371 |
| 70 | controlling | 63.79 | 20.00 | 7.589 |

- From about 36 bar of supply the consumers get 20 bar.
- The flow still rises a little with the supply: a larger let-down cools
  the gas more, and colder gas is denser.
- Below about 36 bar the transmission line cannot deliver the demand at 20
  bar: the regulator is wide open and the consumers get less.

## Limitations

- **No hydrate or ice check here.**  This gas has no water component; for a
  gas with water, the hydrate curve of compositional tutorial 08 applies to
  the regulator outlet.
- **The wide-open valve is a K-factor fitting** (K = 5 unless set); no droop
  or valve curve.
- **Not compared with a published station.**  Isenthalpic expansion is
  checked against the fluid's own h(P, T) in the tests.

## Run

```bash
python src/angelica/tutorials/steady_compositional/11_city_gate_regulator/run.py
```

The GUI file `11_city_gate_regulator.gui.json` holds Part A without the
heater (one thermal segment per pipe, so its outlet temperature differs
from the script's by 0.4 K).
