# Tutorial 01 — Compositional Single Pipe

A single gas pipe transporting a methane/ethane binary mixture from a
high-pressure inlet to a low-pressure outlet.  The EOS flash is performed
at each outer iteration using the `thermo` library, showing how density
changes with pressure while transport properties remain nearly constant
for a gas well above its dew point.

## Network

```
[Node 1] ──── Pipe · 10 km · D = 0.15 m ──── [Node 2]
 P = 100 bar                                   P = 20 bar
 T = 60 °C
```

## Fluid

| Property | Value |
|---|---|
| Components | methane, ethane |
| Composition | CH₄ 80 mol%, C₂H₆ 20 mol% |
| EOS backend | `thermo` library (Peng-Robinson) |
| Two-phase handling | no-slip (homogeneous) |

## EOS flash preview

| Property | Inlet (100 bar) | Outlet (20 bar) |
|---|---|---|
| ρ (kg/m³) | ≈ 68 | ≈ 14 |
| μ (Pa·s × 10⁻⁵) | ≈ 1.18 | ≈ 1.18 |
| Cp (J/(kg·K)) | ≈ 2180 | ≈ 2180 |

Density follows the pressure ratio (~5:1) as expected for a real gas.
Viscosity and Cp are nearly pressure-independent at these conditions.

## Run

```bash
pip install angelica[compositional]   # install thermo if not already done
python tutorials/steady_compositional/01_single_pipe/run.py
```

## Expected output (abridged)

```
EOS flash preview
                            Inlet (100 bar)  Outlet (20 bar)
ρ  (kg/m³)                           68.04           13.61

Converged:         True
Outer iterations:  3
Mass flow:         14.54 kg/s
Velocity @ inlet:  12.09 m/s
Velocity @ outlet: 60.44 m/s   (gas expands 5× as P drops)
Global mass balance error: 0.0000 %
```

Gas velocity nearly quintuples from inlet to outlet because the density
drops 5:1 while mass flow is conserved.
