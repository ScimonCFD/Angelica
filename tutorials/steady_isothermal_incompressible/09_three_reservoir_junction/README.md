# Tutorial 09 — Three-Reservoir Junction

This tutorial implements the classical **three-reservoir junction problem** from
hydraulic engineering: a single junction node supplied by one reservoir and
draining to two others at lower head.

## What This Case Is

The three-reservoir junction problem is a standard benchmark for any pipe-flow
solver. Three reservoirs with fixed hydraulic heads are connected at a common
junction through individual pipes. The unknowns are:

- the head at the junction, H_J
- the flow rate through each of the three pipes

Because the pipes operate in the turbulent regime, the friction factor depends
on Reynolds number (Colebrook-White closure), making the system implicitly
nonlinear. The solver must simultaneously determine H_J and the compatible
friction factors.

The problem has an analytic structure: if H_J is guessed, the three pipe flows
can be computed directly and continuity at J checked. The correct H_J satisfies
mass balance exactly. The Newton-Raphson pressure-correction scheme converges
to this value without an outer guess loop.

## Network

```
Reservoir A (30 m)
      |
      | Pipe A-J: D=0.200m, L=1000m
      |
  Junction J  ─── Pipe J-B: D=0.150m, L=800m ───  Reservoir B (20 m)
      |
      | Pipe J-C: D=0.100m, L=600m
      |
Reservoir C (5 m)
```

| Node | Type     | Head (m) | Pressure (Pa)   |
|------|----------|----------|-----------------|
| 1    | Source   | 30.00    | 293 784.98      |
| 2    | Junction | —        | computed        |
| 3    | Sink     | 20.00    | 195 856.65      |
| 4    | Sink     | 5.00     | 48 964.16       |

| Pipe | Diameter (m) | Length (m) | Roughness (m) |
|------|--------------|------------|---------------|
| A→J  | 0.200        | 1 000      | 0.0001        |
| J→B  | 0.150        | 800        | 0.0001        |
| J→C  | 0.100        | 600        | 0.0001        |

Fluid: water at 998.25 kg/m³, 0.001 Pa·s.

## Expected Results

| Quantity                  | Value           |
|---------------------------|-----------------|
| Junction head H_J         | 25.26 m         |
| Flow A→J                  | 112.35 m³/h     |
| Flow J→B                  | 62.36 m³/h      |
| Flow J→C                  | 49.99 m³/h      |
| Mass balance error at J   | < 1 × 10⁻⁴ kg/s |

The junction head (25.26 m) lies between the heads of reservoirs B (20 m) and
A (30 m), which is required for flow to move from A to J and from J to both B
and C simultaneously.

## What This Tutorial Demonstrates

- Three pressure-fixed boundaries with no pre-specified flow direction
- Junction head determination via pressure-correction Newton-Raphson
- Colebrook-White turbulent friction with implicit friction factor
- Mass conservation across a Y-junction

## Run

```bash
python3 tutorials/steady_isothermal_incompressible/09_three_reservoir_junction/run.py
```

Open in the GUI:

```bash
PYTHONPATH=src python3 -m angelica.gui.app
```

then `File → Open → three_reservoir_junction.gui.json`, press **Run**.
