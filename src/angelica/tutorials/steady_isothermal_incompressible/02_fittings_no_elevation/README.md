# Tutorial 02 — Fittings, No Elevation

Same network as tutorial 01, but with two 90° fittings (K = 1.5) inserted on
two of the pipe connections. No elevation change.

## Network layout

Same 6-node topology as tutorial 01. Fittings are inserted as intermediate
nodes between the pipe segments:

| Device           | Start → End | Diameter (m) | Length or K |
|-----------------|-------------|--------------|-------------|
| Pipe 1a         | 1 → 7       | 0.05         | 150 m       |
| Fitting (K=1.5) | 7 → 8       | 0.05         | —           |
| Pipe 1b         | 8 → 2       | 0.05         | 50 m        |
| Pipe 2a         | 2 → 9       | 0.025        | 10 m        |
| Fitting (K=1.5) | 9 → 10      | 0.025        | —           |
| Pipe 2b         | 10 → 3      | 0.025        | 20 m        |
| Pipe 3          | 2 → 4       | 0.025        | 40 m        |
| Pipe 4          | 4 → 3       | 0.025        | 15 m        |
| Pipe 5          | 5 → 4       | 0.05         | 300 m       |
| Pipe 6          | 3 → 6       | 0.05         | 60 m        |

| Node | Type   | Condition  |
|------|--------|------------|
| 1    | Source | 251 300 Pa |
| 5    | Source | 201 300 Pa |
| 6    | Sink   | 101 300 Pa |

## What this tutorial demonstrates

- Effect of local loss fittings on a multi-loop network
- Fitting K-coefficient model (minor loss)
- Comparison baseline for tutorials 03–06 (same topology)

## Run

```bash
python3 tutorials/steady_isothermal_incompressible/02_fittings_no_elevation/run.py
```
