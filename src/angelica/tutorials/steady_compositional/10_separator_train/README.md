# Tutorial 10 — A separator train (two-stage separation)

A volatile oil with free water flows through a 5 km flowline to a
high-pressure (HP) separator.
- The HP gas goes to export and the water to disposal.
- The oil goes through its level-control valve (LCV) to a low-pressure (LP)
  separator.
- The LP separator's gas goes to a flare and its liquid to the stock tank.

```
Manifold ─ flowline 5 km, 6" ─ HP separator ─┬─ gas line 1 km, 10" ─ Gas export (P_exp)
120 bar, 70 °C                               ├─ water LCV ─ line ─ Disposal 5 bar
                                             └─ oil LCV ─ oil line 500 m, 6" ─ LP separator ─┬─ flare line 100 m, 8" ─ Flare (P_flare)
                                                                                               └─ liquid LCV ─ tank line ─ Stock tank 1.2 bar
```

Feed (mole fractions): C1 0.62, C2 0.09, C3 0.07, nC4 0.05, nC5 0.03,
nC6 0.03, nC10 0.08, water 0.03.

## How it is modelled

- **Separators.** Both are ordinary nodes: their pressures follow from the
  export and flare lines.  (A separator can also be a pressure boundary,
  its gas then leaving through the boundary.)
- **Phase split.** Each separator flashes the arriving streams at its
  pressure and temperature (Peng-Robinson, with free water by the
  immiscible-water model).
- **Stock-tank volumes.** A separator's "oil" rate is the volume of its oil
  phase flashed once to standard conditions (1 atm, 15.56 °C). So:
  - the HP separator's oil rate is what the stock tank would get from the
    HP separator alone (one stage before the tank);
  - the LP separator's oil rate is what it gets after the second stage.
- **Valves.** The liquid crosses each LCV at constant enthalpy, so the oil
  cools as it flashes into the LP separator (46.9 → 37.3 °C).

## Results

**Base case (export 30 bar, flare 5 bar).**

| | P (bar) | T (°C) | Feed (kg/s) | Gas (kg/s) | Oil (kg/s) | Water (kg/s) | Oil (m³/d) |
|---|---|---|---|---|---|---|---|
| HP separator | 31.4 | 46.9 | 42.14 | 19.29 | 22.26 | 0.59 | 2 373.6 |
| LP separator | 5.02 | 37.3 | 22.26 | 1.69 | 20.57 | 0 | 2 491.0 |

- The LP stage raises the stock-tank liquid from 2 374 to 2 491 m³/d.
- Mass and energy balances close (energy residual 0.03 kW).
- LCV drops: water 26.4 bar, HP oil 22.9 bar, LP liquid 3.0 bar.

**HP separator pressure.** Stock-tank liquid per tonne of feed (m³/t):

| P_HP (bar) | One stage | Two stages |
|---|---|---|
| 11.8 | 0.6543 | 0.6647 |
| 18.1 | **0.6589** | 0.6779 |
| 31.4 | 0.6519 | **0.6841** |
| 50.7 | 0.6343 | 0.6781 |
| 80.2 | 0.6068 | 0.6601 |

There is an optimum separator pressure: about 18 bar with the HP separator
alone, about 30 bar with the LP stage.

**LP stage pressure (export 30 bar), set by the flare pressure.** Gain over
one stage:

| P_LP (bar) | T_LP (°C) | Gain |
|---|---|---|
| 2.09 | 29.4 | +2.5 % |
| 5.02 | 37.3 | +4.9 % |
| 10.00 | 42.1 | +4.1 % |

## Limitations

- **Untranslated densities in the flow model.** Peng-Robinson liquid
  densities in the lines are not volume-translated; the standard volumes are
  (see compositional tutorial 09).
- **Not compared with a published separator test.** The split is the PR
  flash.

## Run

```bash
python src/angelica/tutorials/steady_compositional/10_separator_train/run.py
```

The GUI file `10_separator_train.gui.json` holds the base case (one thermal
segment per pipe, so its numbers differ from the script's by up to 0.4 %).
