# Tutorial 15 — A back-pressure regulator holds the separator pressure

The network of tutorial 14 with a back-pressure regulator (BPR) at the start
of the separator's gas line.  The regulator holds the separator at its set
pressure whatever the export pressure downstream, as long as it can.

```
Reservoir 1 ─ Well ─ WH1 ─ flowline 3 km ─┐                    ┌─ BPR (P_set) ─ gas line 4 km, 6" ─ Gas export (P_exp)
                                          ├─ Manifold ─ trunk ─ Separator
Reservoir 2 ─ Well ─ WH2 ─ flowline 2 km ─┘    8 km, 6"        ├─ oil LCV ─ oil line 3 km, 4"  ─ Stock tank 0.2 MPa
                                                               └─ water LCV ─ water line 1 km, 3" ─ Disposal 0.3 MPa
```

## How a regulator is modelled

```python
PressureRegulator(component_id="BPR", start_node=SEPARATOR, end_node=..., diameter_m=0.1524,
                  set_pressure_pa=3.0e6, mode="upstream")
```

- `mode="upstream"` (back-pressure regulator) holds the start node;
  `mode="downstream"` (pressure-reducing regulator) holds the end node
  (compositional tutorial 11).
- The solver puts each regulator in one of three states:
  - **controlling**: the held node is at the set pressure; the regulator
    passes whatever arrives (here, the separator's gas);
  - **open**: the set point cannot be held — here, the export is above it —
    and the regulator is a fitting with its wide-open K (`loss_coefficient`,
    5 by default, an assumption);
  - **closed**: nothing arrives to hold the node up (e.g. the wells dead).
- Flow passes from start to end only.  The fluid crosses it at constant
  enthalpy.  `result.regulators` gives each regulator's state, pressures and
  flow.

## Results

**Part A — export 1.5 MPa.**

| | P_sep (MPa) | Oil (m³/d) | Gas (kSm³/d) | Regulator |
|---|---|---|---|---|
| No regulator (tutorial 14) | 1.537 | 705.5 | 83.4 | — |
| BPR at 3 MPa | 3.000 | 622.3 | 72.6 | controlling, drop 1.47 MPa |

**Part B — export pressure (set point 3 MPa).**

| Export (MPa) | P_sep (MPa) | State | Oil (m³/d) |
|---|---|---|---|
| 0.5 | 3.000 | controlling | 622.3 |
| 1.5 | 3.000 | controlling | 622.3 |
| 2.5 | 3.000 | controlling | 622.3 |
| 3.0 | 3.012 | open | 621.4 |
| 3.5 | 3.509 | open | 565.6 |
| 5.0 | 5.004 | open | 459.5 |

- Below 3 MPa the separator, and so the wells, do not see the export
  pressure at all.
- Above it the regulator is wide open, and the results are tutorial 14's
  (459.5 m³/d at 5 MPa in both).
- At 0.5 MPa export, tutorial 14 does not converge (its oil line falls in a
  Beggs-Brill flow-pattern gap); here the separator stays at 3 MPa and it
  does.

**Part C — set point (export 1 MPa).**

| P_set (MPa) | Oil (m³/d) | Gas (kSm³/d) |
|---|---|---|
| 1.5 | 707.2 | 83.7 |
| 2.0 | 682.5 | 80.3 |
| 3.0 | 622.3 | 72.6 |
| 4.0 | 533.1 | 63.9 |
| 5.0 | 459.9 | 55.4 |

A higher separator pressure costs oil from these wells.  In black-oil it
gains nothing at the stock tank (tutorial 14); the stage effect that can
justify a higher pressure is in compositional tutorial 10.

## Limitations

- **The wide-open valve is a K-factor fitting** (K = 5 unless set).  No
  droop, lock-up or valve curve.
- **Not compared with a published case.**  The states are checked against
  independent solves: controlling against the network on the other side
  solved alone, open against a check valve of the same K.
- **With every flow zero** (both wells dead and the export above the
  separator), the separator is isolated and its pressure is not determined.

## Run

```bash
python src/angelica/tutorials/steady_black_oil/15_back_pressure_regulator/run.py
```

The GUI file `15_back_pressure_regulator.gui.json` holds Part A with the
regulator: the separator's gas connection starts with a "Regulator"
component (holds: upstream, set pressure 3 MPa) followed by its pipe.
