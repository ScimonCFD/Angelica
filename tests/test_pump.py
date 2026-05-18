from __future__ import annotations

import sys
import unittest
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[1] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.closures.pump import PumpCurveModel
from angelica.core.components import Pump
from angelica.core.state import NodeState, PumpState


class PumpCurveModelTests(unittest.TestCase):
    def test_single_point_epanet_curve_hits_design_point(self) -> None:
        model = PumpCurveModel()
        density = 998.25
        design_head_pa = -density * 9.81 * 60.96

        link_state = PumpState(
            component=Pump(
                start_node=1,
                end_node=2,
                diameter_m=0.30,
                curve_points_q_head=((227.124707, 60.96),),
            ),
            start_node=NodeState(node_id=1, pressure_pa=0.0),
            end_node=NodeState(node_id=2, pressure_pa=-design_head_pa),
        )

        velocity = model.calculate_velocity(
            link_state,
            delta_p=design_head_pa,
            density=density,
            viscosity=0.001,
        )
        volumetric_flow_m3_per_h = velocity * link_state.area_m2 * 3600.0

        self.assertAlmostEqual(volumetric_flow_m3_per_h, 227.124707, places=3)

    def test_single_point_epanet_curve_shuts_off_above_shutoff_head(self) -> None:
        model = PumpCurveModel()
        density = 998.25
        shutoff_head_pa = -density * 9.81 * 81.4

        link_state = PumpState(
            component=Pump(
                start_node=1,
                end_node=2,
                diameter_m=0.30,
                curve_points_q_head=((227.124707, 60.96),),
            ),
            start_node=NodeState(node_id=1, pressure_pa=0.0),
            end_node=NodeState(node_id=2, pressure_pa=-shutoff_head_pa),
        )

        velocity = model.calculate_velocity(
            link_state,
            delta_p=shutoff_head_pa,
            density=density,
            viscosity=0.001,
        )
        coupling = model.calculate_coupling(link_state, density=density, viscosity=0.001)

        self.assertEqual(velocity, 0.0)
        self.assertEqual(coupling, 0.0)


if __name__ == "__main__":
    unittest.main()
