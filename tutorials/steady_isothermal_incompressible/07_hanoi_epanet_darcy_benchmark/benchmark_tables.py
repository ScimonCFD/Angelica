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


TERMINAL_NODE_MAP = {
    13: 112,
    22: 121,
}

EXPECTED_LINK_FLOWS_M3_PER_H = {
    1: 19940.00,
    2: 19050.00,
    3: 8011.47,
    4: 7881.47,
    5: 7156.47,
    6: 6151.47,
    7: 4801.47,
    8: 4251.47,
    9: 3726.47,
    10: 2000.00,
    11: 1500.00,
    12: 940.00,
    13: 1201.47,
    14: 586.47,
    15: 306.47,
    16: 242.30,
    17: 1107.30,
    18: 2452.30,
    19: 2512.30,
    20: 7676.23,
    21: 1415.00,
    22: 485.00,
    23: 4986.23,
    24: 3239.90,
    25: 2419.90,
    26: -1031.23,
    27: -131.23,
    28: 238.77,
    29: 701.33,
    30: 411.33,
    31: 51.33,
    32: -308.67,
    33: 413.67,
    34: 1218.67,
}

EXPECTED_NODE_HEADS_M = {
    1: 100.00,
    2: 98.97,
    3: 86.32,
    4: 81.51,
    5: 75.56,
    6: 69.36,
    7: 63.07,
    8: 55.78,
    9: 50.05,
    10: 45.89,
    11: 41.32,
    12: 38.02,
    13: 28.21,
    14: 34.45,
    15: 32.70,
    16: 32.52,
    17: 51.52,
    18: 72.82,
    19: 78.58,
    20: 75.52,
    21: 45.89,
    22: 32.29,
    23: 51.05,
    24: 46.20,
    25: 37.08,
    26: 28.08,
    27: 27.44,
    28: 43.60,
    29: 28.23,
    30: 28.00,
    31: 28.66,
    32: 35.34,
}


def main() -> None:
    scene_path = Path(__file__).with_name("hanoi_epanet_darcy_benchmark.gui.json")
    scene = load_scene_from_file(scene_path)
    case = build_network_case_from_scene(scene)
    solver = build_solver_from_scene(scene)
    result = solver.solve(case)

    rho = case.fluid_model.density_kg_per_m3
    gravity = 9.81

    print(f"Case: {case.name}")
    print(f"Converged: {result.converged}")
    print()
    print("Trunk Node Heads (m)")
    print("node   expected   netSim    delta")
    for benchmark_node_id in range(1, 33):
        result_node_id = TERMINAL_NODE_MAP.get(benchmark_node_id, benchmark_node_id)
        head_m = result.node_pressures_pa[result_node_id] / (rho * gravity)
        expected_head_m = EXPECTED_NODE_HEADS_M[benchmark_node_id]
        delta = head_m - expected_head_m
        print(
            f"node {benchmark_node_id:>2}: "
            f"{expected_head_m:8.3f} {head_m:8.3f} {delta:8.3f}"
        )

    print()
    print("Trunk Link Flows (m^3/h)")
    print("link   expected   netSim    delta")
    for component, flow_result in zip(case.components[:34], result.component_flows[:34]):
        link_id = int(component.component_id.split("_")[1])
        expected_flow_m3_per_h = EXPECTED_LINK_FLOWS_M3_PER_H[link_id]
        calc_flow_m3_per_h = flow_result.volumetric_flow_m3_per_h
        delta = calc_flow_m3_per_h - expected_flow_m3_per_h
        print(
            f"link {link_id:>2}: "
            f"{expected_flow_m3_per_h:8.3f} {calc_flow_m3_per_h:8.3f} {delta:8.3f}"
        )


if __name__ == "__main__":
    main()
