from __future__ import annotations

from angelica.cases import build_steady_water_network_inlet_flow_boundary_case
from angelica.core.settings import SolverSettings
from angelica.io.reporting import print_solve_result
from angelica.solvers import SteadyIsothermalIncompressibleSolver


def main() -> None:
    case = build_steady_water_network_inlet_flow_boundary_case()
    solver = SteadyIsothermalIncompressibleSolver(
        SolverSettings(
            turbulent_iterations=60,
            pressure_relaxation=1.0,
            colebrook_residual_tolerance=1e-4,
            pressure_correction_abs_tolerance_pa=1e-3,
            pressure_correction_rel_tolerance=1e-6,
        )
    )
    result = solver.solve(case)
    print_solve_result(result)


if __name__ == "__main__":
    main()
