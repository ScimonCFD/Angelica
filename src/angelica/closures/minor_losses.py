from __future__ import annotations

import math

from .pressure_drop import PressureDropCorrelation

_K_MAX = 1.0e9  # sentinel for K=∞ presets (e.g. closed check valve)


class MinorLossModel(PressureDropCorrelation):
    def calculate_velocity(
        self,
        link_state,
        delta_p: float,
        density: float,
        viscosity: float,
        tolerance: float | None = None,
    ) -> float:
        K = link_state.component.loss_coefficient
        if math.isinf(K):
            return 0.0
        if delta_p > 0.0:
            return math.sqrt(2.0 * delta_p / (K * density))
        return -math.sqrt(-2.0 * delta_p / (K * density))

    def calculate_coupling(self, link_state, density: float, viscosity: float) -> float:
        K = link_state.component.loss_coefficient
        if math.isinf(K):
            return -link_state.area_m2 / _K_MAX  # near-zero coupling — effectively blocked
        return -2.0 * link_state.area_m2 / (
            K * max(abs(link_state.velocity_m_per_s), 1e-12)
        )

