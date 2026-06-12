# Tutorial 11 — Crude Oil Pipeline

This tutorial demonstrates how to simulate a **crude oil gathering pipeline**
using Angelica's built-in dead-oil property module.  Fluid density and
viscosity are computed automatically from API gravity and operating temperature
via the **Beggs & Robinson (1975)** correlation — no manual look-up required.

## What This Case Shows

- How to use `dead_oil_density_kg_per_m3` and `dead_oil_viscosity_pa_s` from
  `angelica.properties` to supply fluid properties from field data (API gravity,
  operating temperature).
- A branched gathering pipeline with one inlet and two pressure outlets.
- Fully turbulent flow (Re = 19 000 – 31 000) solved with the Colebrook-White
  friction factor.

## Fluid Properties — Beggs & Robinson (1975)

| Property         | Value                                     |
|------------------|-------------------------------------------|
| Crude type       | Dead oil, 32°API                          |
| Temperature      | 65 °C (isothermal — used for μ only)      |
| Density          | 864.6 kg/m³  (from API gravity)           |
| Viscosity        | 4.254 cP  (0.004254 Pa·s)                |

API-gravity-to-density conversion:

    SG = 141.5 / (API + 131.5)
    ρ  = SG × 999.064  kg/m³

Beggs & Robinson dead-oil viscosity:

    log₁₀(log₁₀(μ_dead + 1)) = (3.0324 − 0.02023 API) − 1.163 log₁₀(T_°F)

## Network

```
P_in = 500 000 Pa ─── Pipe 1 (D=0.10 m, L=1 000 m) ─── Node 2 ─┬─ Pipe 2 (D=0.08 m, L=500 m) ─── P_out = 101 325 Pa
  (Node 1)                                                        │
                                                                  └─ Pipe 3 (D=0.06 m, L=300 m) ─── P_out = 101 325 Pa
```

| Node | Type              | Pressure (Pa) |
|------|-------------------|---------------|
| 1    | Pressure inlet    | 500 000       |
| 2    | Interior junction | (solved)      |
| 3    | Pressure outlet A | 101 325       |
| 4    | Pressure outlet B | 101 325       |

| Pipe | From → To | Diameter (m) | Length (m) | Roughness (m)       |
|------|-----------|--------------|------------|---------------------|
| 1    | 1 → 2     | 0.10         | 1 000      | 4.6 × 10⁻⁵ (steel) |
| 2    | 2 → 3     | 0.08         | 500        | 4.6 × 10⁻⁵ (steel) |
| 3    | 2 → 4     | 0.06         | 300        | 4.6 × 10⁻⁵ (steel) |

## Expected Results

| Quantity                     | Value          |
|------------------------------|----------------|
| Pipe 1 flow (m³/h)           | 43.00          |
| Pipe 2 flow (m³/h)           | 26.69          |
| Pipe 3 flow (m³/h)           | 16.31          |
| Node 2 pressure (bar)        | 2.547          |
| Pipe 1 Re                    | ≈ 30 900       |
| Pipe 2 Re                    | ≈ 24 000       |
| Pipe 3 Re                    | ≈ 19 500       |
| Turbulent iterations         | 28             |
| Converged                    | True           |

All pipes operate in the fully turbulent regime (Re ≫ 2 300).  Flow splits
between Pipe 2 and Pipe 3 are governed entirely by pressure-drop balance at
Node 2.  Although Pipe 3 is shorter (300 m vs 500 m), its smaller diameter
(0.06 m vs 0.08 m) dominates: hydraulic resistance scales as L/D⁵, giving
Pipe 3 roughly 2.5× the resistance of Pipe 2 and therefore the lower flow
rate.

## How to Run

```bash
# Quick run — fluid properties + solver output
python3 tutorials/steady_isothermal_incompressible/11_crude_oil_pipeline/run.py
```

Open in the GUI:

```bash
PYTHONPATH=src python3 -m angelica.gui.app
```

then `File → Open → crude_oil_pipeline.gui.json`, press **Run**.

Alternatively, set the fluid interactively:

1. Open the material dialog (toolbar or **Edit → Material**).
2. Select **Crude oil**.
3. Enter **API gravity = 32** and **Temperature = 65 °C**.
4. Click **Apply** — density and viscosity are computed automatically.
