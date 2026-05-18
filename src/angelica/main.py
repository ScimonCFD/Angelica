from __future__ import annotations

from .cases import build_steady_water_network_case
from .core.settings import SolverSettings
from .io.reporting import print_solve_result
from .solvers import SteadyIsothermalIncompressibleSolver


def build_default_solver() -> SteadyIsothermalIncompressibleSolver:
    settings = SolverSettings(
        turbulent_iterations=60,
        pressure_relaxation=1.0,
        colebrook_residual_tolerance=1e-4,
        pressure_correction_abs_tolerance_pa=1e-3,
        pressure_correction_rel_tolerance=1e-8,
    )
    return SteadyIsothermalIncompressibleSolver(settings)


def main() -> None:
    case = build_steady_water_network_case()
    solver = build_default_solver()
    result = solver.solve(case)
    print_solve_result(result)


if __name__ == "__main__":
    main()
