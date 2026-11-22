# Tutorial 12 — Wells with production tubing

Each well carries 2 100 m of production tubing rising 2 000 m from the bottom
hole to the wellhead.  The `Well` component combines the reservoir inflow
(composite Vogel IPR) and the tubing (Beggs-Brill with the local black-oil
PVT), so the network sees each well through its deliverability curve
P_wh(q).

## Network

```
Reservoir 1 ─ Well (IPR + tubing) ─ WH1 ─ flowline 3 km ─┐
Reservoir 2 ─ Well (IPR + tubing) ─ WH2 ─ flowline 2 km ─┼─ Manifold ─ trunk 8 km ─ Separator (1.5 MPa)
Reservoir 3 ─ Well (IPR + tubing) ─ WH3 ─ flowline 4 km ─┘
```

| Well | Crude | P_r (MPa) | T (°C) | J (m³/d/bar) |
|---|---|---|---|---|
| 1 | 40 °API, GOR 150, WOR 0.1 | 25 | 85 | 6 |
| 2 | 32 °API, GOR 40, WOR 0.3 | 22 | 80 | 4 |
| 3 | 25 °API, GOR 15, WOR 0.8 | 17 | 75 | 3 |

## Why a Well component

In a vertical two-phase pipe the liquid holdup depends on the velocity, so
the tubing pressure drop is not monotonic in the flow rate (the J-shaped
tubing performance curve): one pressure drop can match several flow rates.
Solved as ordinary network pipes, the risers make the pressure correction
jump between solutions and the network does not converge, whatever the
relaxation.  The `Well` component integrates the tubing for each rate, keeps
the stable branch of the deliverability curve (the right-hand intersection
of nodal analysis) and hands the network a single-valued, monotonic law
ṁ(P_wh).  If the network asks for a wellhead pressure above the curve's
maximum the well is dead and produces nothing.

## Results (3.0 in tubing)

| Well | P_wh (MPa) | P_wf (MPa) | Oil (m³/d) | |
|---|---|---|---|---|
| 1 | 4.08 | 15.18 | 547 | producing |
| 2 | 2.69 | 17.74 | 170 | producing |
| 3 | 2.62 | — | 0 | dead: column heavier than P_r |

The network converges in 4 outer iterations.  Well 1's deliverability curve
rises from 11.19 MPa (static) to 11.43 MPa at 76 m³/d before falling — gas
lightens the column at low rates.  Larger tubing raises total production:
505 m³/d (50.3 mm), 629 m³/d (62.0 mm), 718 m³/d (76.2 mm).

The tubing is isothermal at the reservoir temperature.  The GUI file
reproduces the 3.0 in case.

## Run

```bash
python src/angelica/tutorials/steady_black_oil/12_wells_with_tubing/run.py
```
