# Tutorial 14 — A three-phase separator in the network

Wells 1 and 2 of tutorial 12 deliver through a trunk line to a production
separator.  The separator is an ordinary node, not a boundary: its pressure
follows from the gas export line downstream.  Its gas, oil and water leave
through their own lines.

```
Reservoir 1 ─ Well ─ WH1 ─ flowline 3 km ─┐                    ┌─ gas line 4 km, 6" ─ Gas export (P_exp)
                                          ├─ Manifold ─ trunk ─ Separator
Reservoir 2 ─ Well ─ WH2 ─ flowline 2 km ─┘    8 km, 6"        ├─ oil LCV ─ oil line 3 km, 4"  ─ Stock tank 0.2 MPa
                                                               └─ water LCV ─ water line 1 km, 3" ─ Disposal 0.3 MPa
```

## How a separator is modelled

A separator is a node with `SeparatorOutlet` links leaving it:

```python
SeparatorOutlet(component_id="oil_LCV", start_node=SEPARATOR, end_node=..., diameter_m=0.0762, phase="oil")
SeparatorOutlet(component_id="water_LCV", start_node=SEPARATOR, end_node=..., diameter_m=0.0508, phase="water")
```

- The streams arriving at the node are mixed and split at its pressure and
  temperature.
  - The oil leaves with the gas still dissolved at the separator, Rs(P, T).
  - The free gas, GOR − Rs, leaves through the node's other links (here the
    gas line) or into its pressure boundary.
  - The water leaves through its own outlet.

  This is the same phase split the black-oil model uses in the pipes.
- `phase` selects the separator type:
  - `"liquid"`: oil + water, a two-phase separator.
  - `"oil"` + `"water"`: three-phase.
  - `"water"` alone: a free-water knockout.
- An outlet link carries its phases whatever its pressure difference: it
  stands for the level-control valve (LCV).
  - Its pressure drop is reported.
  - A negative drop means the liquid could not reach its outlet without a
    pump.
- `result.separators` gives each separator's feed, phase masses and standard
  rates (gas Sm³/d, stock-tank oil Sm³/d, water m³/d), and each outlet's LCV
  pressure drop.

## Results

**Base case (export 1.5 MPa).**
- The separator is at 1.537 MPa and 35.8 °C.
- Of 8.89 kg/s of feed:

  | Phase | Mass (kg/s) | Rate |
  |---|---|---|
  | Gas | 0.82 | 83 400 Sm³/d |
  | Oil | 6.85 | 705.5 Sm³/d |
  | Water | 1.21 | 102 m³/d |

- The stock-tank oil equals the wells' oil.
- LCV drops: 0.50 MPa (oil), 1.23 MPa (water).
- Mass and energy balances close (residual 0.002 kW).

**Export pressure:**

| Export (MPa) | P_sep (MPa) | Oil (m³/d) | Free gas (kSm³/d) | GOR feed | Rs sep. | Well 2 oil (m³/d) |
|---|---|---|---|---|---|---|
| 0.5 | 0.701 | 735.1 | 88.2 | 123 | 3.1 | 179.2 |
| 1.5 | 1.537 | 705.5 | 83.4 | 125 | 7.1 | 158.1 |
| 3.0 | 3.012 | 621.4 | 72.5 | 132 | 15.3 | 101.4 |
| 5.0 | 5.004 | 459.5 | 55.4 | 150 | 29.0 | 1.8 |

- At lower pressure less gas stays dissolved in the oil, and the wells
  produce more.
- At 5 MPa Well 2 (GOR 40) almost dies, so the feed's GOR rises to Well 1's.

**A 2-inch oil line.** The oil would have to leave the LCV at 10.9 MPa,
against a separator at 1.54 MPa. The solution reports that the liquid
needs a pump.

## Limitations

- **No stage-separation effect.** In black-oil the stock-tank oil is fixed
  by the wells' GOR and API whatever the separator pressure.  Compositional
  tutorial 10 shows the stage effect.
- **One gas gravity.** Free and dissolved gas have the same gravity (a
  black-oil assumption).
- **Phantom oil in gas outlets.** A gas outlet carries the gas-only black-oil
  composition (GOR 10⁶), which holds about 0.1 % of its mass as oil.
- **Not compared with a published separator test.** The split is the
  black-oil PVT's own.

## Run

```bash
python src/angelica/tutorials/steady_black_oil/14_separator/run.py
```

The GUI file `14_separator.gui.json` holds the base case.  The separator is
the junction "Separator": each liquid connection starts with a "Separator"
component (the LCV) followed by its pipe.
