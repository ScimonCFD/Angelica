from __future__ import annotations

import sys
import unittest
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[1] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

TUTORIALS_ROOT = SRC_ROOT / "angelica" / "tutorials"

from angelica.gui.io import build_network_case_from_scene, build_solver_from_scene, load_scene_from_file


M_PER_FT = 0.3048
PSI_PER_FT_OF_WATER = 0.433527504
M3H_PER_GPM = 0.22712470704

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


class PumpBenchmarkTests(unittest.TestCase):
    @staticmethod
    def _benchmark_path() -> Path:
        return (
            TUTORIALS_ROOT
            / "steady_isothermal_incompressible"
            / "08_epanet_pump_benchmark"
            / "epanet_pump_benchmark.gui.json"
        )

    def test_epanet_pump_benchmark_matches_reference(self) -> None:
        scene = load_scene_from_file(self._benchmark_path())
        case = build_network_case_from_scene(scene)
        result = build_solver_from_scene(scene).solve(case)

        self.assertTrue(result.converged)

        rho = case.fluid_model.density_kg_per_m3
        gravity = 9.81

        for node_id, expected_pressure_psi in EXPECTED_NODE_PRESSURES_PSI.items():
            pressure_head_ft = (result.node_pressures_pa[node_id] / (rho * gravity)) / M_PER_FT
            pressure_psi = pressure_head_ft * PSI_PER_FT_OF_WATER
            self.assertAlmostEqual(pressure_psi, expected_pressure_psi, delta=0.25)

        for node_id, expected_head_ft in EXPECTED_NODE_HEADS_FT.items():
            pressure_head_ft = (result.node_pressures_pa[node_id] / (rho * gravity)) / M_PER_FT
            total_head_ft = pressure_head_ft + NODE_ELEVATIONS_FT[node_id]
            self.assertAlmostEqual(total_head_ft, expected_head_ft, delta=0.2)

        trunk_flows = {
            int(component.component_id.split("_")[1]): flow_result.volumetric_flow_m3_per_h / M3H_PER_GPM
            for component, flow_result in zip(case.components, result.component_flows)
            if int(component.component_id.split("_")[1]) in EXPECTED_LINK_FLOWS_GPM
        }
        for link_id, expected_flow_gpm in EXPECTED_LINK_FLOWS_GPM.items():
            self.assertAlmostEqual(trunk_flows[link_id], expected_flow_gpm, delta=1.5)


if __name__ == "__main__":
    unittest.main()
