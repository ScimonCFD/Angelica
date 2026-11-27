# Tutorial 16 — A heater station set by its outlet temperature

A 50 km, 8 in. line carries a heavy crude (22 °API, GOR 25, WOR 0.5) that
leaves the field at 80 °C.  A heater station half-way is specified the way
an operator specifies it — by the temperature the oil must leave it at —
and the solver finds the duty.

```
Inlet 10 MPa, 80 °C ─ pipe 25 km ─ [heater: outlet 80 °C, max duty] ─ pipe 25 km ─ Outlet 2 MPa
```

Ground at 10 °C unless stated, U = 3 W/(m²·K), 20 thermal segments per pipe.

## How it is modelled

- `HeatSource(..., outlet_temperature_c=80.0)`: the duty is the enthalpy the
  stream needs to reach 80 °C, ṁ·(h_out − h_in), with the black-oil enthalpy
  model (the same h the energy balance uses).
- `max_duty_w`: beyond it the heater works at its maximum and the oil leaves
  colder than the set point; the result says "duty limited".
- Without `outlet_temperature_c` the heater has a fixed duty (`power_w`).
- `SolveResult.heaters` gives each heater's state, duty and temperatures.

## Results

**Part A — with and without the station** (ambient 10 °C).

| | Flow (kg/s) | T before (°C) | T after (°C) | T outlet (°C) | Duty (MW) |
|---|---|---|---|---|---|
| no heating | 50.45 | 60.90 | 60.92 | 47.35 | 0 |
| station 80 °C | 51.62 | 61.32 | 80.00 | 61.48 | 2.668 |

- The second half runs about 16 K warmer (mean 70 against 54 °C); the
  mixture viscosity there is about 3.5 against 6 mPa·s, and the same 8 MPa
  moves 2.3 % more oil.
- Energy balance residual 0.034 kW of 2 668 kW added (2.9·10⁻⁴ %).

**Part B — ambient temperature**: the duty to hold 80 °C rises from 2.071 MW
at 25 °C to 2.866 MW at 5 °C, about 40 kW per kelvin.

**Part C — a 2.5 MW station**: it holds 80 °C down to 15 °C ambient (2.469
MW); at 10 and 5 °C it is duty limited and the oil leaves at 78.84 and
77.45 °C.

## Limitations

- **Validation.**  The duty is checked against the energy balance and, for
  water, against NIST enthalpies (`tests/test_heater_outlet.py`); the
  black-oil enthalpy model itself is validated only for (∂h/∂P)_T of gas and
  water (see the roadmap).
- **No heater thermal model**: no firing efficiency, no process-side
  approach temperature; the duty is what the stream receives.

## Run

```bash
python src/angelica/tutorials/steady_black_oil/16_heater_station/run.py
```

The GUI file `16_heater_station.gui.json` holds Part A with the station.  Its
pipes use one thermal segment (the analytical exponential profile), so the
duty differs from the script's 20-segment value by 0.7 % (2.687 MW).
