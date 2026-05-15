from __future__ import annotations

from dataclasses import dataclass
import math


FITTING_PRESET_LIBRARY = {
    "regular_90_flanged": {"name": "90 elbow, regular, flanged", "loss_coefficient": 0.3},
    "regular_90_threaded": {"name": "90 elbow, regular, threaded", "loss_coefficient": 1.5},
    "long_radius_90_flanged": {"name": "90 elbow, long radius, flanged", "loss_coefficient": 0.2},
    "long_radius_90_threaded": {"name": "90 elbow, long radius, threaded", "loss_coefficient": 0.7},
    "long_radius_45_flanged": {"name": "45 elbow, long radius, flanged", "loss_coefficient": 0.2},
    "regular_45_threaded": {"name": "45 elbow, regular, threaded", "loss_coefficient": 0.4},
    "return_bend_180_flanged": {"name": "180 return bend, flanged", "loss_coefficient": 0.2},
    "return_bend_180_threaded": {"name": "180 return bend, threaded", "loss_coefficient": 1.5},
    "tee_line_flow_flanged": {"name": "Tee line flow, flanged", "loss_coefficient": 0.2},
    "tee_line_flow_threaded": {"name": "Tee line flow, threaded", "loss_coefficient": 0.9},
    "tee_branch_flow_flanged": {"name": "Tee branch flow, flanged", "loss_coefficient": 1.0},
    "tee_branch_flow_threaded": {"name": "Tee branch flow, threaded", "loss_coefficient": 2.0},
    "union_threaded": {"name": "Union, threaded", "loss_coefficient": 0.08},
    "globe_valve_fully_open": {"name": "Globe valve, fully open", "loss_coefficient": 10.0},
    "angle_valve_fully_open": {"name": "Angle valve, fully open", "loss_coefficient": 2.0},
    "gate_valve_fully_open": {"name": "Gate valve, fully open", "loss_coefficient": 0.15},
    "gate_valve_quarter_closed": {"name": "Gate valve, 1/4 closed", "loss_coefficient": 0.26},
    "gate_valve_half_closed": {"name": "Gate valve, 1/2 closed", "loss_coefficient": 2.1},
    "gate_valve_three_quarters_closed": {
        "name": "Gate valve, 3/4 closed",
        "loss_coefficient": 17.0,
    },
    "swing_check_forward_flow": {"name": "Swing check, forward flow", "loss_coefficient": 2.0},
    "swing_check_backward_flow": {
        "name": "Swing check, backward flow",
        "loss_coefficient": float("inf"),
    },
    "ball_valve_fully_open": {"name": "Ball valve, fully open", "loss_coefficient": 0.05},
    "ball_valve_one_third_closed": {"name": "Ball valve, 1/3 closed", "loss_coefficient": 5.5},
    "ball_valve_two_thirds_closed": {"name": "Ball valve, 2/3 closed", "loss_coefficient": 210.0},
}


@dataclass(frozen=True)
class PressureChanger:
    start_node: int
    end_node: int
    diameter_m: float
    component_id: str = ""

    @property
    def area_m2(self) -> float:
        return math.pi * self.diameter_m**2 / 4.0


@dataclass(frozen=True)
class Pipe(PressureChanger):
    length_m: float = 0.0
    absolute_roughness_m: float = 0.0
    hazen_williams_c: float = 130.0
    height_change_m: float = 0.0


@dataclass(frozen=True)
class Fitting(PressureChanger):
    loss_coefficient: float = 0.0


@dataclass(frozen=True)
class Pump(PressureChanger):
    curve_points_q_head: tuple[tuple[float, float], ...] = ()
