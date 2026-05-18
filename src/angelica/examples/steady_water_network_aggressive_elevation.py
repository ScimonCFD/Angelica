from __future__ import annotations

from angelica.cases import build_steady_water_network_aggressive_elevation_case
from angelica.io.reporting import print_solve_result
from angelica.main import build_default_solver


def main() -> None:
    case = build_steady_water_network_aggressive_elevation_case()
    solver = build_default_solver()
    result = solver.solve(case)
    print_solve_result(result)


if __name__ == "__main__":
    main()
