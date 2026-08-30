from __future__ import annotations

import math

from .gravity import GRAVITY_M_PER_S2
from .pressure_drop import PressureDropCorrelation


class PumpCurveModel(PressureDropCorrelation):
    def calculate_velocity(
        self,
        link_state,
        delta_p: float,
        density: float,
        viscosity: float,
        tolerance: float | None = None,
    ) -> float:
        volumetric_flow_m3_per_h, _slope_q_per_head = self._volumetric_flow_and_slope(
            link_state.component.curve_points_q_head,
            self._pump_head_from_delta_p(delta_p, density),
        )
        volumetric_flow_m3_per_s = volumetric_flow_m3_per_h / 3600.0
        return volumetric_flow_m3_per_s / link_state.area_m2

    def calculate_coupling(self, link_state, density: float, viscosity: float) -> float:
        delta_p = link_state.start_node.pressure_pa - link_state.end_node.pressure_pa
        _flow_q, slope_q_per_head = self._volumetric_flow_and_slope(
            link_state.component.curve_points_q_head,
            self._pump_head_from_delta_p(delta_p, density),
        )
        return slope_q_per_head / (3600.0 * GRAVITY_M_PER_S2)

    @staticmethod
    def _pump_head_from_delta_p(delta_p: float, density: float) -> float:
        return -delta_p / (density * GRAVITY_M_PER_S2)

    def _volumetric_flow_and_slope(
        self,
        raw_points_q_head: tuple[tuple[float, float], ...],
        required_head_m: float,
    ) -> tuple[float, float]:
        if len(raw_points_q_head) == 1:
            return self._single_point_epanet_curve_response(raw_points_q_head, required_head_m)
        return self._piecewise_linear_curve_response(raw_points_q_head, required_head_m)

    def _single_point_epanet_curve_response(
        self,
        single_point_q_head: tuple[tuple[float, float], ...],
        required_head_m: float,
    ) -> tuple[float, float]:
        design_flow_m3_per_h, design_head_m = single_point_q_head[0]
        if design_flow_m3_per_h <= 0.0 or design_head_m <= 0.0:
            raise ValueError("A one-point pump curve requires positive design flow and head.")

        shutoff_head_m = 1.33 * design_head_m
        if required_head_m > shutoff_head_m:
            return 0.0, 0.0

        a_constant = shutoff_head_m
        q_design = design_flow_m3_per_h
        q_max = 2.0 * q_design
        h_design = design_head_m
        h_max = 0.0

        exponent = math.log((a_constant - h_design) / (a_constant - h_max)) / math.log(
            q_design / q_max
        )
        resistance = (a_constant - h_design) / (q_design**exponent)

        head_deficit = max(a_constant - required_head_m, 0.0)
        if head_deficit == 0.0:
            return 0.0, 0.0

        volumetric_flow_m3_per_h = (head_deficit / resistance) ** (1.0 / exponent)
        if volumetric_flow_m3_per_h <= 0.0:
            return 0.0, 0.0

        slope_q_per_head = -volumetric_flow_m3_per_h / (exponent * max(head_deficit, 1e-12))
        return volumetric_flow_m3_per_h, slope_q_per_head

    def _piecewise_linear_curve_response(
        self,
        points_q_head: tuple[tuple[float, float], ...],
        required_head_m: float,
    ) -> tuple[float, float]:
        if len(points_q_head) == 1:
            return points_q_head[0][0], 0.0

        first_q, first_head = points_q_head[0]
        if required_head_m > first_head:
            return 0.0, 0.0

        for index in range(1, len(points_q_head)):
            previous_q, previous_head = points_q_head[index - 1]
            current_q, current_head = points_q_head[index]
            if current_head <= required_head_m <= previous_head:
                slope_q_per_head = (current_q - previous_q) / (current_head - previous_head)
                volumetric_flow_m3_per_h = previous_q + slope_q_per_head * (required_head_m - previous_head)
                return volumetric_flow_m3_per_h, slope_q_per_head

        previous_q, previous_head = points_q_head[-2]
        current_q, current_head = points_q_head[-1]
        slope_q_per_head = (current_q - previous_q) / (current_head - previous_head)
        volumetric_flow_m3_per_h = current_q + slope_q_per_head * (required_head_m - current_head)
        return max(volumetric_flow_m3_per_h, 0.0), slope_q_per_head
