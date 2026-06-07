from __future__ import annotations

from angelica.cases.three_reservoir_junction import build_three_reservoir_junction_case
from angelica.closures import ColebrookPipeCorrelation
from angelica.core.settings import SolverSettings
from angelica.io.reporting import print_solve_result
from angelica.solvers import SteadyIsothermalIncompressibleSolver


def build_example_case():
    return build_three_reservoir_junction_case()


def main() -> None:
    case = build_example_case()
    solver = SteadyIsothermalIncompressibleSolver(
        settings=SolverSettings(
            turbulent_iterations=300,
            pressure_relaxation=0.5,
            pressure_correction_abs_tolerance_pa=1e-3,
            pressure_correction_rel_tolerance=1e-8,
        ),
        turbulent_pipe_correlation=ColebrookPipeCorrelation(),
    )
    result = solver.solve(case)
    print_solve_result(result)


if __name__ == "__main__":
    main()
