from __future__ import annotations

import math

from .pressure_drop import PressureDropCorrelation

# K_FLOOR: minimum loss coefficient for a "transparent" device (pressure_drop_pa = 0).
# Chosen so ΔP ≈ 250 Pa at 7 m/s water flow — negligible vs typical bar-scale systems.
_K_FLOOR = 0.01

# Maximum coupling magnitude (kg/s/Pa).  Caps the heater coupling when v→0 so the
# pressure matrix stays well-conditioned despite the small loss coefficient.
# At 1 kg/s/Pa the condition number vs a typical pipe coupling (≈1e-4) is ~1e4,
# giving floating-point errors well below the 1e-3 Pa convergence tolerance.
_MAX_COUPLING_ABS = 1.0


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
        raw = -2.0 * link_state.area_m2 / (K * max(abs(link_state.velocity_m_per_s), 1e-12))
        # Clamp to avoid ill-conditioning when v→0 with small K
        return max(raw, -_MAX_COUPLING_ABS)
