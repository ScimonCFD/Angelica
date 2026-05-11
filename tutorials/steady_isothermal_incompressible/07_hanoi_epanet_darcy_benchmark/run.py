from __future__ import annotations

import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[3] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from netSim.gui.io import (
    build_network_case_from_scene,
    build_solver_from_scene,
    load_scene_from_file,
)


def main() -> None:
    scene_path = Path(__file__).with_name("hanoi_epanet_darcy_benchmark.gui.json")
    scene = load_scene_from_file(scene_path)
    case = build_network_case_from_scene(scene)
    solver = build_solver_from_scene(scene)
    result = solver.solve(case)

    print(f"Case: {case.name}")
    print(f"Converged: {result.converged}")
    print(f"Laminar iterations: {len(result.laminar_metrics)}")
    print(f"Turbulent iterations: {len(result.turbulent_metrics)}")

    if result.turbulent_metrics:
        final_metrics = result.turbulent_metrics[-1]
        print(
            "Final metrics: "
            f"abs_dp={final_metrics.pressure_correction_abs_pa:.6e} Pa, "
            f"rel_dp={final_metrics.pressure_correction_rel:.6e}, "
            f"max_mass_imbalance={final_metrics.max_nodal_mass_imbalance_kg_per_s:.6e} kg/s"
        )


if __name__ == "__main__":
    main()
