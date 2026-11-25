from __future__ import annotations

import math

from .pressure_drop import PressureDropCorrelation

_K_MAX = 1.0e9  # sentinel for K=∞ presets (e.g. closed check valve)


def check_valve_closed(link_state) -> bool:
    """True for a CheckValve whose downstream pressure is not below upstream."""
    if not getattr(link_state.component, "blocks_reverse_flow", False):
        return False
    p_start = link_state.start_node.pressure_pa
    p_end = link_state.end_node.pressure_pa
    if p_start is None or p_end is None:
        return False
    return float(p_start) <= float(p_end)


def _reject_well_inflow(link_state) -> None:
    comp = link_state.component
    if getattr(comp, "is_well_inflow", False) or getattr(comp, "is_well_model", False):
        raise ValueError(
            "WellInflow / Well links need a black-oil or compositional fluid: solve the "
            "case with SteadyBlackOilSolver or SteadyCompositionalSolver."
        )
    if getattr(comp, "is_separator_outlet", False):
        raise ValueError(
            "SeparatorOutlet links need a black-oil or compositional fluid: solve the "
            "case with SteadyBlackOilSolver or SteadyCompositionalSolver."
        )
    if getattr(comp, "is_pressure_regulator", False):
        raise ValueError(
            "PressureRegulator links need SteadyBlackOilSolver or SteadyCompositionalSolver "
            "(a single-phase gas or liquid can be given as a compositional fluid)."
        )


def k_fitting_velocity(link_state, delta_p: float, density: float) -> float:
    """Velocity of a K-factor fitting: ΔP = K·ρ·v²/2 (zero when blocked)."""
    K = link_state.component.loss_coefficient
    if math.isinf(K) or check_valve_closed(link_state):
        return 0.0
    if delta_p > 0.0:
        return math.sqrt(2.0 * delta_p / (K * density))
    return -math.sqrt(-2.0 * delta_p / (K * density))


def k_fitting_coupling(link_state) -> float:
    """Secant dṁ/dΔP of a K-factor fitting (near zero when blocked)."""
    K = link_state.component.loss_coefficient
    if math.isinf(K) or check_valve_closed(link_state):
        return -link_state.area_m2 / _K_MAX  # near-zero coupling — effectively blocked
    return -2.0 * link_state.area_m2 / (
        K * max(abs(link_state.velocity_m_per_s), 1e-12)
    )


class MinorLossModel(PressureDropCorrelation):
    def calculate_velocity(
        self,
        link_state,
        delta_p: float,
        density: float,
        viscosity: float,
        tolerance: float | None = None,
    ) -> float:
        _reject_well_inflow(link_state)
        return k_fitting_velocity(link_state, delta_p, density)

    def calculate_coupling(self, link_state, density: float, viscosity: float) -> float:
        _reject_well_inflow(link_state)
        return k_fitting_coupling(link_state)

