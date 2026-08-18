# Tutorial 03 — Fittings with Elevation

Same network as tutorial 02 (fittings, K = 1.5) but with elevation changes
applied to every pipe segment. This tests the hydrostatic pressure term.

## Network layout

Identical topology to tutorial 02. Elevation changes added per pipe:

| Pipe segment  | Height change (m) |
|---------------|-------------------|
| 1 → 7         | −2.0 (downhill)   |
| 8 → 2         | −3.0 (downhill)   |
| 2 → 9         | +4.0 (uphill)     |
| 10 → 3        | +3.0 (uphill)     |
| 2 → 4         | +2.0 (uphill)     |
| 4 → 3         | +5.0 (uphill)     |
| 5 → 4         | −3.0 (downhill)   |
| 3 → 6         | +4.0 (uphill)     |

Boundary conditions same as tutorial 02:

| Node | Type   | Condition  |
|------|--------|------------|
| 1    | Source | 251 300 Pa |
| 5    | Source | 201 300 Pa |
| 6    | Sink   | 101 300 Pa |

## What this tutorial demonstrates

- Hydrostatic pressure contribution (ρ g Δz) for pipes with elevation change
- Combined friction + minor loss + gravity pressure drop

## Run

```bash
python3 tutorials/steady_isothermal_incompressible/03_fittings_with_elevation/run.py
```
