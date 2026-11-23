# Tutorial 09 — Gas and gas-condensate wells

A lean-gas well and a gas-condensate well produce through their tubing,
flowlines and a common trunk line to a separator at 60 bar.  The `Well`
component (IPR + tubing, black-oil tutorial 12) works here with the
Peng-Robinson equation of state.

## Network

```
Reservoir G ─ gas well (3 000 m)        ─ WH-G ─ flowline 3 km, 6" ─┐
                                                                    ├─ Manifold ─ trunk 10 km, 8" ─ Separator (60 bar)
Reservoir C ─ condensate well (3 200 m) ─ WH-C ─ flowline 4 km, 6" ─┘
```

| Well | Fluid | P_r | T_r | IPR (gas back-pressure) |
|---|---|---|---|---|
| Gas | 86 % methane, 0.5 % n-decane | 250 bar | 100 °C | C = 150, n = 0.80 |
| Condensate | 74 % methane, 4 % n-decane | 280 bar | 110 °C | C = 120, n = 0.75 |

Both: 3-in. tubing, Gray correlation, Ramey heat loss to the formation.

## How a compositional well is computed

- **IPR.** The back-pressure law gives the gas rate in standard m³/d.  It
  becomes a mass rate through a flash of the stream to standard conditions
  (101 325 Pa, 15.56 °C): kg of stream per standard m³ of gas.  The same
  flash gives the condensate (stock-tank liquid) volume, with the Peneloux
  volume correction.  The oil-rate IPRs (linear, Vogel, Fetkovich) work for
  oil and condensate streams; a dry gas has no stock-tank liquid and needs
  the back-pressure law.
- **Tubing.** Each well's stream is flashed once on a pressure–temperature
  table (about 900 points, 5–10 s) and interpolated along the tubing.
  Against a march that flashes at every step the table gives the wellhead
  pressure within 0.05 % and temperature within 0.01 K.
- **Temperature.** Isenthalpic inflow, then the first law along the tubing
  with Ramey's heat loss.  The gas cools on the drawdown (Joule-Thomson),
  unlike the black-oil liquids of tutorial 12, which warm slightly.

## Results (separator 60 bar)

| Well | P_wh | T_wh | Gas | Condensate |
|---|---|---|---|---|
| Gas | 70.1 bar | 44.1 °C | 649 kSm³/d | 22 Sm³/d |
| Condensate | 69.7 bar | 59.7 °C | 434 kSm³/d | 169 Sm³/d |

Energy balance of the network: residual below 0.001 kW.

Along the condensate well's tubing the liquid mass fraction grows from 26 %
to 34 % (condensation as the fluid climbs and cools) while the liquid volume
fraction falls from 9 % to 5 % (the gas expands faster).

### Separator pressure (Part D)

| Separator | Gas well | Condensate well |
|---|---|---|
| 60 bar | 6.40 kg/s | 6.23 kg/s |
| 150 bar | 3.76 kg/s | 3.47 kg/s |
| 186 bar | 1.32 kg/s | dead |

The condensate well dies first although its reservoir pressure is higher:
its liquid-laden column is heavier (shut-in wellhead pressure 175 bar
against the gas well's 190 bar), and its curve peaks at 184.5 bar.

## Limitations

- The flow model uses Peng-Robinson densities without volume translation
  (the Peneloux correction is applied only to standard volumes and surface
  tension).  For this condensate's liquid in the tubing the difference is
  below 1 % (70–150 bar); for heavy liquids it is larger — 12 % for the
  n-C16-rich stock-tank oil of the volatile-oil comparison in the ROADMAP
  log — and the column comes out lighter than with translated densities.
- Not validated against a published gas-well example.

## Run

```bash
python src/angelica/tutorials/steady_compositional/09_gas_condensate_wells/run.py
```
