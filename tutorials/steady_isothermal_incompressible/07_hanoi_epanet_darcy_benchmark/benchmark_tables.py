from __future__ import annotations

import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[3] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.gui.io import (
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
    3: 8011.36,
    4: 7881.36,
    5: 7156.36,
    6: 6151.36,
    7: 4801.36,
    8: 4251.36,
    9: 3726.36,
    10: 2000.00,
    11: 1500.00,
    12: 940.00,
    13: 1201.36,
    14: 586.36,
    15: 306.36,
    16: 243.22,
    17: -1108.22,
    18: -2453.22,
    19: -2513.22,
    20: 7675.42,
    21: 1415.00,
    22: 485.00,
    23: 4985.42,
    24: 3238.91,
    25: 2418.91,
    26: -1030.42,
    27: -130.42,
    28: 239.58,
    29: 701.51,
    30: 411.51,
    31: 51.51,
    32: -308.49,
    33: 413.49,
    34: 1218.49,
}

EXPECTED_NODE_HEADS_M = {
    1: 100.00,
    2: 99.04,
    3: 87.07,
    4: 82.31,
    5: 76.42,
    6: 70.20,
    7: 64.28,
    8: 57.22,
    9: 51.58,
    10: 47.42,
    11: 42.79,
    12: 39.36,
    13: 29.13,
    14: 36.50,
    15: 34.69,
    16: 34.49,
    17: 53.58,
    18: 74.16,
    19: 79.84,
    20: 76.33,
    21: 48.60,
    22: 36.04,
    23: 52.71,
    24: 47.78,
    25: 38.78,
    26: 30.05,
    27: 29.38,
    28: 45.15,
    29: 29.91,
    30: 29.65,
    31: 30.32,
    32: 36.93,
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
    print("node   expected   Angelica    delta")
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
    print("link   expected   Angelica    delta")
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
