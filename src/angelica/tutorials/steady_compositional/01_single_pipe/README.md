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

Flash at 60 °C.  Reference values in brackets: CoolProp (GERG-2008 type
Helmholtz mixture model with its reference viscosity correlations).

| Property | Inlet (100 bar) | Outlet (20 bar) |
|---|---|---|
| ρ (kg/m³) | 81.9 (80.1) | 14.2 (14.1) |
| μ (Pa·s × 10⁻⁵) | 1.46 (1.48) | 1.22 (1.22) |
| Cp (J/(kg·K)) | 2869 (2880) | 2294 (2290) |

Density follows the pressure ratio (~5:1) as expected for a real gas.  Gas
viscosity rises by ~20 % between 20 and 100 bar: Angelica adds a dense-gas
(Jossi-Stiel-Thodos) term to the low-pressure mixture viscosity.  Cp also
rises with pressure, as the gas departs from ideal behaviour.

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
