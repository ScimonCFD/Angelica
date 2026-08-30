# Angelica

Open-source steady-state pipe network simulator for isothermal, non-isothermal, compressible, black-oil, and compositional flows.

## Installation

```bash
pip install angelica
```

For the graphical interface:

```bash
pip install "angelica[gui]"
```

## Quick Start

```python
import angelica as ag

# Define network topology
case = ag.NetworkCase(
    name="simple_pipe",
    components=(
        ag.Pipe(start_node=1, end_node=2, diameter_m=0.1, length_m=1000.0,
                roughness_m=4.6e-5),
    ),
    pressure_inlets=(ag.PressureBoundary(node_id=1, pressure_pa=6e5),),
    pressure_outlets=(ag.PressureBoundary(node_id=2, pressure_pa=1e5),),
    fluid=ag.SingleComponentFluid(density_kg_per_m3=998.0, viscosity_pa_s=1e-3),
)

result = ag.SteadyIsothermalIncompressibleSolver().solve(case)
ag.print_solve_result(result)
```

## Simulation Modes

| Mode | Solver | Fluid |
|------|--------|-------|
| Isothermal incompressible | `SteadyIsothermalIncompressibleSolver` | `SingleComponentFluid` |
| Non-isothermal incompressible | `SteadyNonIsothermalIncompressibleSolver` | `ThermalFluid` |
| Compressible (gas) | `SteadyCompressibleSolver` | `CompressibleFluid` |
| Black-oil | `SteadyBlackOilSolver` | `BlackOilFluid` |
| Compositional | `SteadyCompositionalSolver` | `CompositionalFluid` |

## Network Components

- **`Pipe`** — straight pipe segment with friction, elevation, and optional heat exchange.
- **`Fitting`** — minor loss element (K-factor or equivalent length).
- **`Pump`** — centrifugal pump defined by a Q–H curve.
- **`HeatSource`** — fixed power input element for thermal networks.
- **`PressureChanger`** — base class for all link components.

## Links

- [Source code](https://github.com/ScimonCFD/Angelica)
- [Issue tracker](https://github.com/ScimonCFD/Angelica/issues)
