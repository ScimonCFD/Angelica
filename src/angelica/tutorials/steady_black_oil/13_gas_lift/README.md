# Tutorial 13 — Gas lift

Well 3 of tutorial 12 is a heavy, wet crude (25 °API, GOR 15, WOR 0.8): its
static column outweighs the reservoir pressure and the well is dead.  Gas
lift injects gas into the tubing through a valve; above the valve the
column is lighter and the well flows.

```
Reservoir ─ Well (IPR + 2 100 m tubing, gas injected at depth) ─ WH ─ flowline 2 km, 4" ─ Separator
17 MPa, 75 °C                                                                              1.5 MPa
```

## How gas lift is modelled

`Well(..., gas_lift_rate_sm3_per_day=60_000, gas_lift_depth_m=1800)`:

- Below the valve the tubing carries the reservoir fluid; above it the
  reservoir fluid plus the lift gas.  The lift gas reaches the valve at the
  formation temperature of that depth and stays gas (no mass transfer with
  the oil; for black-oil, valid when the oil at the valve is saturated, as
  here).
- The lift gas leaves the well at the wellhead node, where the network
  receives it as an internal source: the flowline carries production + lift
  gas, with a gas-oil ratio that includes the lift gas.  The link's flow and
  `result.well_states[...].mass_flow_kg_per_s` are the reservoir's
  production; `gas_lift_kg_per_s` is the injected gas.
- Lift gas gravity (black-oil) or mole fractions (compositional) are
  optional; by default the produced gas.

## Results

Gas-lift performance curve (valve at 1 800 m):

| Lift gas (Sm³/d) | Oil (m³/d) | P_wh (MPa) | Valve P (MPa) |
|---|---|---|---|
| 0 | 0 (dead) | — | — |
| 10 000 | 9.9 | 1.54 | 14.10 |
| 20 000 | 120.9 | 1.78 | 10.39 |
| 40 000 | 146.9 | 2.04 | 9.52 |
| 70 000 | 156.7 | 2.37 | 9.19 |
| 100 000 | 155.5 | 2.65 | 9.23 |
| 250 000 | 117.4 | 3.82 | 10.51 |

The oil rate peaks near 70 000 Sm³/d: beyond it the extra gas adds more
friction (in the tubing and the flowline) than it removes weight.  Below
about 15 000 Sm³/d the well barely flows.

Injection depth at 60 000 Sm³/d: 37 m³/d at 600 m, 108 at 1 200 m, 155 at
1 800 m, 174 at the bottom (2 100 m) — a deeper valve lightens more of the
column but needs more injection pressure at the valve (2.8 → 11.2 MPa).

Energy balance of the network: residual 0.002 kW, with the lift gas
entering at the valve and its potential energy from the valve to the
wellhead as explicit terms.

## Not validated

Against published gas-lift performance data (none accessible).  The
annulus is not modelled: the surface injection pressure is not computed.

## Run

```bash
python src/angelica/tutorials/steady_black_oil/13_gas_lift/run.py
```
