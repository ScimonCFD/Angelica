# Solvers

Angelica provides five steady-state solvers. All implement `BaseSolver.solve(case)` and return a `SolveResult`.

## Isothermal Incompressible

```python
from angelica import SteadyIsothermalIncompressibleSolver, SolverSettings

solver = SteadyIsothermalIncompressibleSolver(
    settings=SolverSettings(max_iterations=200, tolerance=1e-6)
)
result = solver.solve(case)
```

Pressure-velocity coupling uses the SIMPLE algorithm with separate laminar and turbulent iteration loops.

## Non-Isothermal Incompressible

```python
from angelica import SteadyNonIsothermalIncompressibleSolver, NonIsothermalSolverSettings

solver = SteadyNonIsothermalIncompressibleSolver(
    non_isothermal_settings=NonIsothermalSolverSettings(
        max_temperature_iterations=50,
        temperature_tolerance=1e-4,
    )
)
result = solver.solve(case)
```

Requires at least one `ThermalBoundary` with `bc_type="fixed_temperature"`.

## Compressible

```python
from angelica import SteadyCompressibleSolver, CompressibleSolverSettings

solver = SteadyCompressibleSolver(
    compressible_settings=CompressibleSolverSettings(max_density_iterations=30)
)
result = solver.solve(case)
```

Density is updated from the EOS at each outer iteration.

## Black-Oil

```python
from angelica import SteadyBlackOilSolver, BlackOilSolverSettings

solver = SteadyBlackOilSolver(
    black_oil_settings=BlackOilSolverSettings(max_pvt_iterations=30)
)
result = solver.solve(case)
```

## Compositional

```python
from angelica import SteadyCompositionalSolver, CompositionalSolverSettings

solver = SteadyCompositionalSolver(
    compositional_settings=CompositionalSolverSettings(max_flash_iterations=30)
)
result = solver.solve(case)
```

Flash calculations use the `thermo` library (Michelsen TP-flash).
