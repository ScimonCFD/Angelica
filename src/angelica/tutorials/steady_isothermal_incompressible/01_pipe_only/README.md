# Tutorial 01 — Pipe-Only Network

A small looped network solved with Darcy-Weisbach / Colebrook-White using
only pipes (no fittings, no elevation change).

## Network layout

```
     [1]──P1──[2]──P2──[3]──P6──[6]
      │         │         │
    251.3 kPa  P3         P4
                │         │
               [4]──P5──[?]
                │
              201.3 kPa [5]
```

More precisely:

| Node | Type   | Condition         |
|------|--------|-------------------|
| 1    | Source | 251 300 Pa        |
| 5    | Source | 201 300 Pa        |
| 6    | Sink   | 101 300 Pa        |
| 2, 3, 4 | Junction | —            |

| Pipe | Start → End | Diameter (m) | Length (m) |
|------|-------------|--------------|------------|
| P1   | 1 → 2       | 0.05         | 200        |
| P2   | 2 → 3       | 0.025        | 30         |
| P3   | 2 → 4       | 0.025        | 40         |
| P4   | 4 → 3       | 0.025        | 15         |
| P5   | 5 → 4       | 0.05         | 300        |
| P6   | 3 → 6       | 0.05         | 60         |

Fluid: water at 998.25 kg/m³, 0.001 Pa·s.
Roughness: 0.000045 m for all pipes.

## What this tutorial demonstrates

- Multi-loop pressure-driven flow with pure Darcy-Weisbach closure
- No additional loss sources (no fittings, no elevation, no pumps)
- Two pressure inlets and one pressure outlet

## Run

```bash
python3 tutorials/steady_isothermal_incompressible/01_pipe_only/run.py
```
