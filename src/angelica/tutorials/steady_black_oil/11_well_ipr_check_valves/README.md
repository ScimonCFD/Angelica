# Tutorial 11 — Wells with IPR and check valves

Three wells deliver oil through their inflow performance relationship (IPR)
into a gathering network with Beggs-Brill flowlines.  A check valve at each
wellhead stops back-flow from the manifold into the weakest well.

## Network

```
Reservoir 1 ─IPR─ WH1 ─check─ flowline 4 km, 4"  ─┐
Reservoir 2 ─IPR─ WH2 ─check─ flowline 2 km, 4"  ─┼─ Manifold ─ trunk 8 km, 6" ─ Separator (1.5 MPa)
Reservoir 3 ─IPR─ WH3 ─check─ flowline 6 km, 3"  ─┘
```

| Well | Crude | P_r (MPa) | T (°C) | J (m³/d/bar) |
|---|---|---|---|---|
| 1 | 40 °API, GOR 150, WOR 0.1 | 8.0 | 80 | 25 |
| 2 | 32 °API, GOR 40, WOR 0.3 | 6.5 | 65 | 20 |
| 3 | 25 °API, GOR 15, WOR 0.8 | 3.2 | 50 | 15 |

Reservoir pressures are shut-in pressures referred to the wellhead: the
tubing is folded into each well's deliverability curve, as is common in
gathering-network models.  Tutorial 12 adds the production tubing with the
`Well` component.

## Components

- `WellInflow` — linear or composite-Vogel IPR between a reservoir node
  (pressure boundary) and the bottom-hole / wellhead node.  The stock-tank
  oil rate is converted to a mass rate with the produced gas and water.
- `CheckValve` — K-factor fitting (K = 2 open) that closes when the
  downstream pressure is higher.

## Results

| | Well 1 | Well 2 | Well 3 | Manifold |
|---|---|---|---|---|
| A — no check valves | 7.32 kg/s | 6.38 kg/s | −0.15 kg/s (injected) | 3.29 MPa |
| B — check valves | 7.31 kg/s | 6.35 kg/s | shut in | 3.31 MPa |

With a linear IPR the same network produces about 8 % more oil than with
Vogel, which accounts for the gas liberated near the wellbore.

The GUI file reproduces case B.

## Run

```bash
python src/angelica/tutorials/steady_black_oil/11_well_ipr_check_valves/run.py
```
