"""
Well inflow performance relationships (IPR) for ``WellInflow`` links.

The stock-tank oil rate q_o follows a linear productivity index or the
composite Vogel (1968) relationship; the mass rate adds the gas and water
produced with each stock-tank cubic metre of oil:

    ṁ = q_o · (ρ_o,sc + GOR·ρ_g,sc + WOR·ρ_w,sc)

References
----------
Vogel, J.V. (1968). Inflow Performance Relationships for Solution-Gas Drive
Wells. JPT 20(1), 83–92.
Brown, K.E. (1984). The Technology of Artificial Lift Methods, Vol. 4 —
composite IPR for undersaturated reservoirs.
"""
from __future__ import annotations

from dataclasses import dataclass

_PA_PER_BAR = 1.0e5
_S_PER_DAY = 86_400.0
_AIR_DENSITY_SC = 1.2250      # kg/m³, as in angelica.properties.black_oil
_WATER_DENSITY_SC = 1_025.0   # kg/m³, as in angelica.properties.black_oil


@dataclass(frozen=True)
class WellInflowPoint:
    reservoir_pressure_pa: float
    flowing_pressure_pa: float
    bubble_point_pa: float          # capped at the reservoir pressure (Vogel only)
    oil_rate_sc_m3_per_day: float   # > 0 production, < 0 injection
    mass_flow_kg_per_s: float       # signed like the oil rate
    mass_per_oil_sc_m3: float       # kg produced per stock-tank m³ of oil


def oil_rate_sc_m3_per_day(
    p_res: float, p_wf: float, pi_m3_per_day_per_bar: float, model: str, p_bubble: float,
) -> float:
    """Stock-tank oil rate from a linear or composite-Vogel IPR (pressures in Pa)."""
    J = pi_m3_per_day_per_bar / _PA_PER_BAR          # m³/d per Pa
    if p_wf >= p_res or model == "linear":
        return J * (p_res - p_wf)
    pb = min(p_bubble, p_res)
    if p_wf >= pb:
        return J * (p_res - p_wf)
    x = p_wf / pb
    return J * (p_res - pb) + J * pb / 1.8 * (1.0 - 0.2 * x - 0.8 * x * x)


def oil_rate_derivative(
    p_res: float, p_wf: float, pi_m3_per_day_per_bar: float, model: str, p_bubble: float,
) -> float:
    """−dq_o/dP_wf (m³/d per Pa), used for the pressure–flow coupling."""
    J = pi_m3_per_day_per_bar / _PA_PER_BAR
    if p_wf >= p_res or model == "linear":
        return J
    pb = min(p_bubble, p_res)
    if p_wf >= pb:
        return J
    x = p_wf / pb
    return J / 1.8 * (0.2 + 1.6 * x)


def mass_per_oil_sc_m3(fluid) -> float:
    """Mass produced (oil + gas + water) per stock-tank m³ of oil."""
    from angelica.properties.dead_oil import dead_oil_density_kg_per_m3

    rho_o = dead_oil_density_kg_per_m3(fluid.api_gravity)
    rho_g = fluid.gas_gravity * _AIR_DENSITY_SC
    return rho_o + fluid.gor_sc_m3_per_m3 * rho_g + fluid.wor_sc_m3_per_m3 * _WATER_DENSITY_SC


def well_inflow_point(link_state, fluid) -> WellInflowPoint:
    """Evaluate a WellInflow link at the current node pressures."""
    comp = link_state.component
    p_res = float(link_state.start_node.pressure_pa)
    p_wf = float(link_state.end_node.pressure_pa)
    T = link_state.start_node.temperature_c
    if T is None:
        T = getattr(fluid, "reference_temperature_c", 60.0)
    p_b = fluid.pvt(p_res, T).bubble_point_pa if comp.ipr_model == "vogel" else p_res
    q = oil_rate_sc_m3_per_day(p_res, p_wf, comp.productivity_index_m3_per_day_per_bar, comp.ipr_model, p_b)
    m_per = mass_per_oil_sc_m3(fluid)
    return WellInflowPoint(
        reservoir_pressure_pa=p_res,
        flowing_pressure_pa=p_wf,
        bubble_point_pa=min(p_b, p_res),
        oil_rate_sc_m3_per_day=q,
        mass_flow_kg_per_s=q * m_per / _S_PER_DAY,
        mass_per_oil_sc_m3=m_per,
    )
