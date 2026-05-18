from __future__ import annotations

from angelica.cases.steady_water_network_no_fittings import build_steady_water_network_no_fittings_case
from angelica.core.settings import SolverSettings
from angelica.io.reporting import print_solve_result
from angelica.solvers import SteadyIsothermalIncompressibleSolver


def build_example_case():
    return build_steady_water_network_no_fittings_case()


def main() -> None:
    case = build_example_case()
    solver = SteadyIsothermalIncompressibleSolver(
        SolverSettings(
            turbulent_iterations=60,
            pressure_relaxation=1.0,
            colebrook_residual_tolerance=1e-4,
            pressure_correction_abs_tolerance_pa=1e-3,
            pressure_correction_rel_tolerance=1e-8,
        )
    )
    result = solver.solve(case)
    print_solve_result(result)


if __name__ == "__main__":
    main()
