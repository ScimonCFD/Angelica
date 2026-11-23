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

| Well | P_wh (MPa) | P_wf (MPa) | Oil (m³/d) | T_wh (°C) | |
|---|---|---|---|---|---|
| 1 | 4.04 | 15.17 | 548 | 70.4 | producing |
| 2 | 2.68 | 18.00 | 160 | 52.7 | producing |
| 3 | 2.61 | — | 0 | — | dead: column heavier than P_r |

Well 1's deliverability curve rises from 10.97 MPa (static) to 11.04 MPa at
76 m³/d before falling — gas lightens the column at low rates.  Larger
tubing raises total production: 498 m³/d (50.3 mm), 620 m³/d (62.0 mm),
708 m³/d (76.2 mm).

### Tubing correlation (Part D)

| Correlation | Well 1 max P_wh | Total oil |
|---|---|---|
| Beggs-Brill | 11.04 MPa | 708 m³/d |
| Hagedorn-Brown (modified) | 11.30 MPa | 769 m³/d |
| Gray | 11.96 MPa | 759 m³/d |

Select it with `Well(..., tubing_correlation="hagedorn_brown")` or in the
GUI well dialog. Hagedorn-Brown is meant for oil wells, Gray for gas and
gas-condensate wells (its holdup does not tend to 1 at zero rate). The
spread between correlations is part of the uncertainty of any well model.

### Tubing temperature (Part E)

The fluid enters the tubing at its reservoir enthalpy and loses heat to the
formation on the way up (Ramey's model; the rock is at the reservoir
temperature at the bottom and 0.03 K/m cooler towards the surface).  Well 1
at 6.7 kg/s:

| Height above bottom (m) | Pressure (MPa) | Fluid (°C) | Formation (°C) |
|---|---|---|---|
| 0 | 15.17 | 89.3 | 85.0 |
| 800 | 10.31 | 84.6 | 61.0 |
| 1600 | 6.06 | 76.0 | 37.0 |
| 2000 | 4.06 | 70.4 | 25.0 |

It gives 195 kW to the formation.  The same well would arrive at 42 °C at
1 kg/s and 64 °C at 4 kg/s: a slow well has time to cool.  The 4 K rise at
the bottom hole is the liquid expanding at constant enthalpy through the
drawdown.  The network's energy balance closes to 0.002 kW with the wells'
formation heat and lifting work as explicit terms.

Thermal inputs of a `Well` (all optional): `geothermal_gradient_k_per_m`,
`formation_conductivity_w_per_mk` (0 = adiabatic tubing),
`formation_diffusivity_m2_per_s`, `production_time_days`,
`wellbore_diameter_m`, `tubing_u_w_per_m2k`.  Not validated against
measured flowing temperatures.

> **The default `tubing_u_w_per_m2k = inf` is a bound, not a best estimate.**
> It means the tubing wall, annulus, casing and cement offer no resistance
> to heat flow: only the rock limits it.  The heat loss is therefore an
> upper bound and the wellhead temperatures above are lower bounds — the
> real well arrives warmer.  Give the overall coefficient of the completion
> (referred to the tubing inner diameter) when it is known; a gas- or
> liquid-filled annulus typically lowers the heat loss substantially.

### Back-pressure (Part F)

Raising the separator pressure lifts the manifold:

| Separator | Manifold | Well 1 oil | Well 2 oil |
|---|---|---|---|
| 1.5 MPa | 2.61 MPa | 548 m³/d | 160 m³/d |
| 4.0 MPa | 4.38 MPa | 493 m³/d | 40 m³/d |
| 6.0 MPa | 6.18 MPa | 415 m³/d | 0 |

At 6 MPa the manifold is above Well 2's maximum wellhead pressure
(5.29 MPa): Well 2 dies against the manifold, its flowline stagnates, and
Well 1 carries all the production.

A well whose curve peaks at a finite rate — a gassy well like Well 1, whose
maximum is at 76 m³/d — cannot flow steadily below that rate.  Between its
static (shut-in) wellhead pressure and the peak it is **bistable**: dead, or
producing on the stable branch.  The solver prefers producing when that is
consistent.  If producing would lift the manifold above the peak, the well
is dead: shut in, it cannot start against the manifold, and it would need
unloading or gas lift.  `result.well_states` gives each well's state and
the reason:

| State | Meaning |
|---|---|
| `producing` | on the stable branch |
| `dead_column` | the static column outweighs the reservoir |
| `dead_back_pressure` | wellhead above the curve's maximum |
| `dead_minimum_rate` | bistable band; producing is not consistent, dead is |
| `unstable` | neither state is consistent (`converged=False`) |

`unstable` needs a dead well whose wellhead falls *below* its static
pressure (so it would start) while its minimum stable rate lifts the
wellhead above the maximum — e.g. a weak well on its own thin flowline.  The
real well then produces intermittently.  (v1.23.1 wrongly said that any well
below its minimum stable rate had no steady state.)

The GUI file reproduces the 3.0 in case.

## Run

```bash
python src/angelica/tutorials/steady_black_oil/12_wells_with_tubing/run.py
```
