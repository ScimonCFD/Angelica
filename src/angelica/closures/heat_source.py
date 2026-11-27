from __future__ import annotations

import math

from .pressure_drop import PressureDropCorrelation

# K_FLOOR: minimum loss coefficient for a "transparent" device (pressure_drop_pa = 0).
# Chosen so ΔP ≈ 250 Pa at 7 m/s water flow — negligible vs typical bar-scale systems.
_K_FLOOR = 0.01

# Velocity floor of the coupling (m/s): keeps the coupling of a nearly
# transparent device finite when v → 0 (at most 2A/(K·1e-4)).  Before
# v1.29.0 the coupling itself was capped at 1 in magnitude; a transparent
# heater in a water line at 1.4 m/s needs 6.4, so the pressure correction
# under-estimated its conductance six-fold and the solve did not converge
# (246 kg/s through the heater of a 64 kg/s line).  A floor of 0.01 m/s
# fixed that but under-estimated the coupling between 1e-4 and 1e-2 m/s,
# and a heater after a dead well then carried a phantom 0.025 kg/s and the
# solve did not converge (found in review); 1e-4 m/s converges both.
_V_FLOOR = 1e-4


class HeatSourceModel(PressureDropCorrelation):
    """Hydraulic model for inline heater/cooler components.

    Rated mode:  ΔP = ΔP_rated · (ṁ/ṁ_rated)² — equivalent to a fixed K factor
                 K = 2·ΔP_rated·ρ·A² / ṁ_rated²

    Fixed mode:  ΔP ≈ pressure_drop_pa (lagged coefficient updated each iteration)
                 K_iter = 2·ΔP_fixed / (ρ·v_prev²)

    Both modes degrade gracefully to a transparent device when pressure_drop_pa = 0.
    """

    def _effective_k(self, link_state, density: float) -> float:
        hs = link_state.component
        dp = hs.pressure_drop_pa
        if dp <= 0.0:
            return _K_FLOOR
        A = link_state.area_m2
        if hs.pressure_drop_mode == "rated":
            mdot_r = hs.rated_mass_flow_kg_per_s
            if mdot_r <= 0.0:
                return _K_FLOOR
            return max(2.0 * dp * density * A ** 2 / mdot_r ** 2, _K_FLOOR)
        # fixed (lagged): derive K from previous iteration velocity
        v = abs(link_state.velocity_m_per_s)
        if v < 1e-12:
            return _K_FLOOR
        return max(2.0 * dp / (density * v ** 2), _K_FLOOR)

    def calculate_velocity(
        self,
        link_state,
        delta_p: float,
        density: float,
        viscosity: float,
        tolerance: float | None = None,
    ) -> float:
        K = self._effective_k(link_state, density)
        if delta_p > 0.0:
            return math.sqrt(2.0 * delta_p / (K * density))
        return -math.sqrt(-2.0 * delta_p / (K * density))

    def calculate_coupling(self, link_state, density: float, viscosity: float) -> float:
        K = self._effective_k(link_state, density)
        return -2.0 * link_state.area_m2 / (K * max(abs(link_state.velocity_m_per_s), _V_FLOOR))
