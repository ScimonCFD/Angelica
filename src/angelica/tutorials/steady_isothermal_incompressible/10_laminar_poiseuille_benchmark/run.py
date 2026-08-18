from __future__ import annotations

import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[3] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.cases.laminar_parallel_pipes import build_laminar_parallel_pipes_case
from angelica.solvers import SteadyIsothermalIncompressibleSolver


def main() -> None:
    case = build_laminar_parallel_pipes_case()
    solver = SteadyIsothermalIncompressibleSolver()
    result = solver.solve(case)

    print(f"Case: {case.name}")
    print(f"Converged: {result.converged}")
    print(f"Laminar iterations: {len(result.laminar_metrics)}")
    print(f"Turbulent iterations: {len(result.turbulent_metrics)}")
    print()
    print("Component Flows")
    print(f"{'Component':<20}  {'Mass flow (kg/s)':>16}  {'Vol. flow (m³/h)':>18}")
    for comp in result.component_flows:
        print(
            f"{comp.label:<20}  "
            f"{comp.mass_flow_kg_per_s:>16.6f}  "
            f"{comp.volumetric_flow_m3_per_h:>18.6f}"
        )

    if result.laminar_metrics:
        final_metrics = result.laminar_metrics[-1]
        print()
        print(
            "Final laminar metrics: "
            f"abs_dp={final_metrics.pressure_correction_abs_pa:.6e} Pa, "
            f"rel_dp={final_metrics.pressure_correction_rel:.6e}, "
            f"max_mass_imbalance_rel={final_metrics.max_nodal_mass_imbalance_rel:.6e}"
        )


if __name__ == "__main__":
    main()
