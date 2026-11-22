# Tutorial 08 — Hydrate risk in a subsea wet-gas line

A wet natural gas with produced water leaves a platform at 90 bar and 40 °C
and flows 30 km along the seabed (4 °C) to shore at 60 bar. Gas hydrates
form wherever the fluid is colder than the hydrate formation temperature at
the local pressure and free water is present. The same line is solved bare
(U = 15 W/m²K) and insulated (U = 1.5 W/m²K).

Fluid (mol %): methane 81, ethane 6, propane 3, n-butane 1, water 9.

## Hydrate model

van der Waals–Platteeuw with the Munck et al. (1988) Langmuir constants,
Peng-Robinson fugacities, liquid-water branch (T ≥ 0 °C), pure water (no
inhibitor). Validated against the measured quadruple points of methane,
ethane and propane hydrates (within 0.8 K). Every compositional pipe reports
`hydrate_temperature_c` and `hydrate_subcooling_k = T_hydrate − T_fluid`
(positive: inside the hydrate region) at its colder end. The GUI draws the
hydrate curve on the phase envelope and labels the links at risk.

## Results

Hydrate curve: 6.9 bar at 0 °C, 31.7 bar at 12 °C, 119 bar at 20 °C
(structure II — the 3 % propane lowers the formation pressure to about a
quarter of pure methane's).

| | Bare line | Insulated line |
|---|---|---|
| Arrival temperature | 1.1 °C | 19.5 °C |
| Hydrate temperature at arrival | 16.4 °C | 16.4 °C |
| First section at risk | km 5–10 | none |

The bare line cools to the seabed temperature within about 10 km — and
below it towards shore, because the expanding gas also cools by
Joule-Thomson effect — and runs up to 15 K inside the hydrate region. The
insulated line arrives 3 K above the hydrate temperature.

The GUI file reproduces the bare line.

## Run

```bash
python src/angelica/tutorials/steady_compositional/08_hydrate_risk_subsea/run.py
```
