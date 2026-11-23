"""
Well inflow performance relationships (IPR) for ``WellInflow`` and ``Well``.

Models (pressures in bar inside the formulas, rates in standard m³/d):

    linear        q_o = J·(P_r − P_wf)
    vogel         composite Vogel: linear above the bubble point P_b, then
                  q_o = J·(P_r − P_b) + J·P_b/1.8·(1 − 0.2·x − 0.8·x²), x = P_wf/P_b
    fetkovich     q_o = C·(P_r² − P_wf²)ⁿ                      (oil wells)
    backpressure  q_g = C·(P_r² − P_wf²)ⁿ                      (gas wells)

J in stock-tank m³ of oil per day per bar; C in standard m³/d per bar^(2n)
(of oil for Fetkovich, of gas for the back-pressure equation); n between 0.5
(fully turbulent) and 1 (laminar, Darcy flow).  With n = 1 Fetkovich's
equation is the pseudo-steady-state solution for a pressure-dependent
productivity proportional to pressure.  When P_wf > P_r the well takes
injection by the same law with the sign reversed.

The mass rate adds what is produced with each standard m³ of the reference
phase.  Black-oil: ṁ = q_o·(ρ_o,sc + GOR·ρ_g,sc + WOR·ρ_w,sc), and for the
back-pressure equation q_o = q_g / GOR (GOR = total produced gas).

References
----------
Vogel, J.V. (1968). Inflow Performance Relationships for Solution-Gas Drive
Wells. JPT 20(1), 83–92.
Brown, K.E. (1984). The Technology of Artificial Lift Methods, Vol. 4 —
composite IPR for undersaturated reservoirs.
Fetkovich, M.J. (1973). The Isochronal Testing of Oil Wells. SPE 4529.
Rawlins, E.L., Schellhardt, M.A. (1936). Back-Pressure Data on Natural Gas
Wells and Their Application to Production Practices. U.S. Bureau of Mines
Monograph 7.
"""
from __future__ import annotations

import math
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


IPR_MODELS = ("linear", "vogel", "fetkovich", "backpressure")


@dataclass(frozen=True)
class IPR:
    """Inflow law of one well.  ``rate_basis`` is the phase whose standard
    volume rate the law gives ("oil" or "gas")."""

    model: str
    productivity_index: float = float("nan")   # J, linear / vogel
    coefficient: float = float("nan")          # C, fetkovich / backpressure
    exponent: float = 1.0                      # n

    @property
    def rate_basis(self) -> str:
        return "gas" if self.model == "backpressure" else "oil"

    @property
    def uses_bubble_point(self) -> bool:
        return self.model == "vogel"


def ipr_of(comp) -> IPR:
    """The IPR of a ``WellInflow`` or ``Well`` component."""
    return IPR(comp.ipr_model, float(comp.productivity_index_m3_per_day_per_bar),
               float(getattr(comp, "ipr_coefficient", float("nan"))),
               float(getattr(comp, "ipr_exponent", 1.0)))


def validate_ipr_fields(owner: str, model: str, J: float, C: float, n: float) -> None:
    """Raise ValueError for an incomplete or inconsistent IPR specification."""
    if model not in IPR_MODELS:
        raise ValueError(f"{owner} ipr_model must be one of {IPR_MODELS} (got {model!r})")
    if model in ("linear", "vogel"):
        if math.isnan(J) or J <= 0.0:
            raise ValueError(f"{owner} with ipr_model={model!r} requires a positive "
                             "productivity_index_m3_per_day_per_bar.")
    else:
        if math.isnan(C) or C <= 0.0:
            raise ValueError(f"{owner} with ipr_model={model!r} requires a positive ipr_coefficient "
                             "(standard m³/d per bar^(2n)).")
        if not 0.5 <= n <= 1.0:
            raise ValueError(f"{owner} ipr_exponent must be between 0.5 and 1 (got {n}).")


def _bar(p: float) -> float:
    return p / _PA_PER_BAR


def ipr_rate(ipr: IPR, p_res: float, p_wf: float, p_bubble: float = float("inf")) -> float:
    """Standard-volume rate of the reference phase (m³/d), > 0 production."""
    if ipr.model in ("linear", "vogel"):
        return oil_rate_sc_m3_per_day(p_res, p_wf, ipr.productivity_index, ipr.model, p_bubble)
    d = _bar(p_res) ** 2 - _bar(p_wf) ** 2
    return math.copysign(ipr.coefficient * abs(d) ** ipr.exponent, d)


def ipr_rate_derivative(ipr: IPR, p_res: float, p_wf: float, p_bubble: float = float("inf")) -> float:
    """−dq/dP_wf (m³/d per Pa).  For n < 1 the tangent is infinite at zero
    drawdown; it is evaluated at a drawdown of at least 0.1 % of P_r there."""
    if ipr.model in ("linear", "vogel"):
        return oil_rate_derivative(p_res, p_wf, ipr.productivity_index, ipr.model, p_bubble)
    pr, pw = _bar(p_res), _bar(p_wf)
    d = abs(pr * pr - pw * pw)
    d_min = pr * pr - (0.999 * pr) ** 2
    d = max(d, d_min)
    pw_eff = math.sqrt(max(pr * pr - d, 0.0)) if pw <= pr else math.sqrt(pr * pr + d)
    return ipr.coefficient * ipr.exponent * d ** (ipr.exponent - 1.0) * 2.0 * pw_eff / _PA_PER_BAR


def ipr_open_flow(ipr: IPR, p_res: float, p_bubble: float = float("inf")) -> float:
    """Absolute open flow (rate at P_wf = 0), m³/d of the reference phase."""
    return ipr_rate(ipr, p_res, 0.0, p_bubble)


def ipr_flowing_pressure(ipr: IPR, q: float, p_res: float, p_bubble: float = float("inf")) -> float:
    """Flowing bottom-hole pressure (Pa) for a production rate q ≥ 0 of the
    reference phase (inverse of the IPR); −1 beyond the absolute open flow."""
    if ipr.model == "linear":
        return p_res - q / (ipr.productivity_index / _PA_PER_BAR)
    if ipr.model == "vogel":
        J_pa = ipr.productivity_index / _PA_PER_BAR
        pb = min(p_bubble, p_res)
        q_b = J_pa * (p_res - pb)
        if q <= q_b:
            return p_res - q / J_pa
        c = 1.0 - (q - q_b) * 1.8 / (J_pa * pb)
        if c < 0.0:
            return -1.0
        return (-0.2 + math.sqrt(0.04 + 3.2 * c)) / 1.6 * pb
    pr2 = _bar(p_res) ** 2
    d = (max(q, 0.0) / ipr.coefficient) ** (1.0 / ipr.exponent)
    if d > pr2:
        return -1.0
    return math.sqrt(pr2 - d) * _PA_PER_BAR


def mass_per_rate_unit(fluid, ipr: IPR) -> float:
    """Mass produced (kg) per standard m³ of the IPR's reference phase."""
    m_oil = mass_per_oil_sc_m3(fluid)
    if ipr.rate_basis == "oil":
        return m_oil
    gor = float(fluid.gor_sc_m3_per_m3)
    if gor <= 0.0:
        raise ValueError("The back-pressure IPR gives a gas rate: it needs a black-oil fluid with GOR > 0.")
    return m_oil / gor


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
    ipr = ipr_of(comp)
    p_b = fluid.pvt(p_res, T).bubble_point_pa if ipr.uses_bubble_point else p_res
    q_ref = ipr_rate(ipr, p_res, p_wf, p_b)
    m_ref = mass_per_rate_unit(fluid, ipr)
    m_oil = mass_per_oil_sc_m3(fluid)
    m_dot = q_ref * m_ref / _S_PER_DAY
    return WellInflowPoint(
        reservoir_pressure_pa=p_res,
        flowing_pressure_pa=p_wf,
        bubble_point_pa=min(p_b, p_res),
        oil_rate_sc_m3_per_day=m_dot * _S_PER_DAY / m_oil,
        mass_flow_kg_per_s=m_dot,
        mass_per_oil_sc_m3=m_oil,
    )


def well_inflow_coupling(link_state, fluid) -> float:
    """dṁ/dP_wf (< 0) of a WellInflow link at the current node pressures."""
    comp = link_state.component
    p_res = float(link_state.start_node.pressure_pa)
    p_wf = float(link_state.end_node.pressure_pa)
    T = link_state.start_node.temperature_c
    if T is None:
        T = getattr(fluid, "reference_temperature_c", 60.0)
    ipr = ipr_of(comp)
    p_b = fluid.pvt(p_res, T).bubble_point_pa if ipr.uses_bubble_point else p_res
    dq = ipr_rate_derivative(ipr, p_res, p_wf, p_b)
    return -max(dq * mass_per_rate_unit(fluid, ipr) / _S_PER_DAY, 1e-12)
