# Tutorial 10 — Laminar Poiseuille Benchmark

This tutorial verifies the **laminar initialisation phase** of the solver against
the analytical Hagen-Poiseuille solution for viscous pipe flow.

## What This Case Tests

Three parallel smooth pipes carry a viscous oil between fixed pressure boundaries.
Because the fluid is highly viscous and the pipes are narrow, all flows are well
within the laminar regime (Re < 2300).  The Angelica laminar solver uses the
Poiseuille closure

    v = D² ΔP / (32 μ L)

which should reproduce the analytical result exactly (to floating-point
precision) once the pressure field has converged.

The solver runs with default settings.  Because all Reynolds numbers are well
below 2300, `ColebrookPipeCorrelation` automatically uses f = 64/Re (the
laminar Poiseuille friction factor) for any turbulent-phase iteration that
touches these pipes, so the result is consistent with the exact closed-form
solution regardless of which phase solves it.

## Network

```
P_in = 50 000 Pa  ──  Pipe 1 (D=0.020m, L=10m)  ──┐
                  ──  Pipe 2 (D=0.015m, L= 8m)  ──┤── P_out = 0 Pa
                  ──  Pipe 3 (D=0.025m, L=12m)  ──┘
```

| Node | Type   | Pressure (Pa) |
|------|--------|---------------|
| 1    | Source | 50 000        |
| 2    | Sink   | 0             |

| Pipe | Diameter (m) | Length (m) | Roughness (m) |
|------|--------------|------------|---------------|
| 1    | 0.020        | 10         | 0.0 (smooth)  |
| 2    | 0.015        | 8          | 0.0 (smooth)  |
| 3    | 0.025        | 12         | 0.0 (smooth)  |

**Fluid:** Oil — density = 880 kg/m³, viscosity = 0.1 Pa·s

## Analytical Solution

The Hagen-Poiseuille equation for fully-developed laminar flow in a circular
pipe gives:

    Q = π D⁴ ΔP / (128 μ L)

With ΔP = 50 000 Pa:

| Pipe | Q (m³/s)       | Q (m³/h)  | Re  |
|------|----------------|-----------|-----|
| 1    | 1.9635 × 10⁻⁴  | 0.7069    | 110 |
| 2    | 7.7632 × 10⁻⁵  | 0.2795    | 58  |
| 3    | 3.9936 × 10⁻⁴  | 1.4377    | 179 |

## Reynolds Numbers — Laminar Regime Confirmation

All three pipes have Re < 2300, confirming the laminar-flow assumption.
Maximum Re ≈ 179 (Pipe 3), far below the turbulent transition threshold.

## Expected Results

| Quantity                         | Value              |
|----------------------------------|--------------------|
| Pipe 1 flow (m³/h)               | 0.70686            |
| Pipe 2 flow (m³/h)               | 0.27948            |
| Pipe 3 flow (m³/h)               | 1.43769            |
| Max relative error vs analytical | < 0.001 %          |
| Converged                        | True               |

## How to Run

```bash
# Quick run — shows solver output
python3 tutorials/steady_isothermal_incompressible/10_laminar_poiseuille_benchmark/run.py

# Benchmark — comparison table vs analytical solution
python3 tutorials/steady_isothermal_incompressible/10_laminar_poiseuille_benchmark/benchmark_tables.py
```

Open in the GUI:

```bash
PYTHONPATH=src python3 -m angelica.gui.app
```

then `File → Open → laminar_poiseuille_benchmark.gui.json`, press **Run**.
