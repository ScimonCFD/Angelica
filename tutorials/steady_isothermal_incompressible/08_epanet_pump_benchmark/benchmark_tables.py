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


M_PER_FT = 0.3048
M3H_PER_GPM = 0.22712470704
PSI_PER_FT_OF_WATER = 0.433527504

NODE_ELEVATIONS_FT = {
    2: 0.0,
    3: 710.0,
    4: 700.0,
    5: 695.0,
    6: 700.0,
    7: 850.0,
}

EXPECTED_NODE_HEADS_FT = {
    2: 893.19,
    3: 879.67,
    4: 874.36,
    5: 872.62,
    6: 872.65,
    7: 855.00,
}

EXPECTED_NODE_PRESSURES_PSI = {
    2: 387.02,
    3: 73.52,
    4: 75.55,
    5: 76.96,
    6: 74.81,
    7: 2.17,
}

EXPECTED_LINK_FLOWS_GPM = {
    1: 1049.81,
    2: 559.25,
    3: 165.56,
    4: 90.56,
    5: -9.44,
    6: 474.81,
    7: 1049.81,
}


def main() -> None:
    scene_path = Path(__file__).with_name("epanet_pump_benchmark.gui.json")
    scene = load_scene_from_file(scene_path)
    case = build_network_case_from_scene(scene)
    solver = build_solver_from_scene(scene)
    result = solver.solve(case)

    rho = case.fluid_model.density_kg_per_m3
    gravity = 9.81

    print(f"Case: {case.name}")
    print(f"Converged: {result.converged}")
    print()
    print("Node Pressures (psi)")
    print("node   expected   netSim    delta")
    for node_id in sorted(EXPECTED_NODE_PRESSURES_PSI):
        pressure_head_m = result.node_pressures_pa[node_id] / (rho * gravity)
        pressure_psi = (pressure_head_m / M_PER_FT) * PSI_PER_FT_OF_WATER
        expected_psi = EXPECTED_NODE_PRESSURES_PSI[node_id]
        delta_psi = pressure_psi - expected_psi
        print(f"node {node_id:>2}: {expected_psi:8.3f} {pressure_psi:8.3f} {delta_psi:8.3f}")

    print()
    print("Node Total Heads (ft)")
    print("node   expected   netSim    delta")
    for node_id in sorted(EXPECTED_NODE_HEADS_FT):
        pressure_head_ft = (result.node_pressures_pa[node_id] / (rho * gravity)) / M_PER_FT
        total_head_ft = pressure_head_ft + NODE_ELEVATIONS_FT[node_id]
        expected_head_ft = EXPECTED_NODE_HEADS_FT[node_id]
        delta_head_ft = total_head_ft - expected_head_ft
        print(f"node {node_id:>2}: {expected_head_ft:8.3f} {total_head_ft:8.3f} {delta_head_ft:8.3f}")

    print()
    print("Trunk Link Flows (gpm)")
    print("link   expected   netSim    delta")
    trunk_results = {int(component.component_id.split('_')[1]): flow_result for component, flow_result in zip(case.components, result.component_flows) if int(component.component_id.split('_')[1]) in EXPECTED_LINK_FLOWS_GPM}
    for link_id in sorted(EXPECTED_LINK_FLOWS_GPM):
        calc_flow_gpm = trunk_results[link_id].volumetric_flow_m3_per_h / M3H_PER_GPM
        expected_flow_gpm = EXPECTED_LINK_FLOWS_GPM[link_id]
        delta_flow_gpm = calc_flow_gpm - expected_flow_gpm
        print(f"link {link_id:>2}: {expected_flow_gpm:8.3f} {calc_flow_gpm:8.3f} {delta_flow_gpm:8.3f}")


if __name__ == "__main__":
    main()
