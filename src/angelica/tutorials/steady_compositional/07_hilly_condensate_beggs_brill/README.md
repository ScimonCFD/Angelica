# Tutorial 07 — Gas-condensate line over hilly terrain: homogeneous vs Beggs-Brill

A wet gas with condensate and produced water flows 20 km over hilly terrain,
between 80 bar at the inlet and 55 bar at the outlet. The same network is
solved with the compositional default (homogeneous no-slip flow) and with
Beggs-Brill, which lets the gas slip past the liquid.

## Network

```
Inlet ── A: 5 km, +80 m ── B: 5 km, −60 m ── C: 5 km, +120 m ── D: 5 km, −40 m ── Outlet
80 bar, 45 °C                                                                   55 bar
```

10 in (254 mm) line, U = 4 W/m²K to soil at 10 °C.

Fluid (mol %): methane 78, ethane 8, propane 5, n-butane 3, n-pentane 2,
n-hexane 2, water 2; Peng-Robinson. The liquid is the condensate plus free
water; water vapour travels with the gas.

## Using Beggs-Brill with a compositional fluid

```python
from angelica.closures import BeggsBrillCorrelation
from angelica.solvers import SteadyCompositionalSolver

result = SteadyCompositionalSolver(
    turbulent_pipe_correlation=BeggsBrillCorrelation()
).solve(case)
```

In the GUI choose *Physics → Define Pressure-Drop Model → Beggs-Brill
(two-phase)*. Each pipe's flash gives the gas and liquid densities,
viscosities and the gas–condensate surface tension (Macleod-Sugden); their no-slip mixture
density equals the homogeneous model's density, so the two models differ
only in slip.

## Results

The section table comes from `result.component_flows[i].multiphase`
(`MultiphaseFlowResult`), also exported to CSV/Excel and shown in the GUI
(right-click a link → *Multiphase Flow*).


| Model | Flow (kg/s) |
|---|---|
| Homogeneous | 30.17 |
| Beggs-Brill | 25.54 (15.3 % less) |

| Section | dz (m) | Pattern | C_L | H_L | ΔP elevation (bar) | ΔP friction (bar) |
|---|---|---|---|---|---|---|
| A | +80 | segregated | 0.020 | 0.140 | 1.14 | 4.08 |
| B | −60 | transition | 0.023 | 0.116 | −0.77 | 5.86 |
| C | +120 | transition | 0.024 | 0.132 | 1.56 | 6.27 |
| D | −40 | transition | 0.022 | 0.106 | −0.44 | 7.30 |

At these low mixture velocities the liquid holdup is five to seven times the
no-slip liquid fraction: the line holds about 125 m³ of liquid, against 23 m³
in the homogeneous model. Uphill sections hold the most; downhill the
liquid drains faster. The extra liquid raises the two-phase friction factor,
which takes most of the 25 bar.

The GUI file reproduces the Beggs-Brill case.

## Run

```bash
python src/angelica/tutorials/steady_compositional/07_hilly_condensate_beggs_brill/run.py
```
