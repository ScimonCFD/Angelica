from __future__ import annotations

import math

from .pressure_drop import PressureDropCorrelation


class MinorLossModel(PressureDropCorrelation):
    def calculate_velocity(
        self,
        link_state,
        delta_p: float,
        density: float,
        viscosity: float,
        tolerance: float | None = None,
    ) -> float:
        loss_coefficient = link_state.component.loss_coefficient
        if delta_p > 0.0:
            return math.sqrt(2.0 * delta_p / (loss_coefficient * density))
        return -math.sqrt(-2.0 * delta_p / (loss_coefficient * density))

    def calculate_coupling(self, link_state, density: float, viscosity: float) -> float:
        return -2.0 * link_state.area_m2 / (
            link_state.component.loss_coefficient
            * max(abs(link_state.velocity_m_per_s), 1e-12)
        )

    def velocity_from_pressure_drop(self, delta_p: float, density: float, loss_coefficient: float) -> float:
        if delta_p > 0.0:
            return math.sqrt(2.0 * delta_p / (loss_coefficient * density))
        return -math.sqrt(-2.0 * delta_p / (loss_coefficient * density))

    def coupling_coefficient(self, area_m2: float, loss_coefficient: float, velocity_m_per_s: float) -> float:
        return -2.0 * area_m2 / (loss_coefficient * max(abs(velocity_m_per_s), 1e-12))
