"""
Well deliverability (IPR + production tubing) for ``Well`` links.

For an oil rate q the flowing bottom-hole pressure follows the IPR; the
tubing pressure drop is integrated from the bottom hole to the wellhead with
Beggs-Brill and the local black-oil PVT.  This gives the deliverability curve
P_wh(q), tabulated once per well and fluid.  The network uses its stable
branch — where P_wh falls as q rises, the right-hand intersection of nodal
analysis — so that the link law ṁ(P_wh) is single-valued and monotonic.
Above the curve's maximum wellhead pressure the well is dead (ṁ = 0).

A vertical two-phase riser solved directly as a network pipe is not
single-valued (the J-shaped tubing performance curve) and makes the SIMPLE
pressure correction jump between solutions; folding the tubing into the well
model avoids that.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache

from .beggs_brill import beggs_brill_gradient, phase_properties_from_pvt, surface_tension_from_pvt
from .well_inflow import mass_per_oil_sc_m3

_PA_PER_BAR = 1.0e5
_S_PER_DAY = 86_400.0
_P_MIN_PA = 1.2e5          # the tubing integration stops below ~1.2 bar
_N_RATES = 60              # points of the tabulated deliverability curve


@dataclass(frozen=True)
class WellCurve:
    """Tabulated deliverability curve of one well."""

    oil_rates_sc_m3_per_day: tuple[float, ...]   # ascending, starting at 0
    wellhead_pressures_pa: tuple[float, ...]
    mass_per_oil_sc_m3: float
    peak_index: int                               # index of the maximum P_wh

    @property
    def max_wellhead_pressure_pa(self) -> float:
        return self.wellhead_pressures_pa[self.peak_index]

    def mass_flow(self, p_wh: float) -> tuple[float, float]:
        """(ṁ, dṁ/dP_wh) on the stable branch at wellhead pressure ``p_wh``."""
        q, p = self.oil_rates_sc_m3_per_day, self.wellhead_pressures_pa
        k = self.mass_per_oil_sc_m3 / _S_PER_DAY
        i0 = self.peak_index
        if p_wh >= p[i0]:
            return 0.0, 0.0                       # dead well
        for i in range(i0, len(q) - 1):
            if p[i + 1] <= p_wh <= p[i]:
                dq_dp = (q[i + 1] - q[i]) / (p[i + 1] - p[i])
                return k * (q[i] + dq_dp * (p_wh - p[i])), k * dq_dp
        # Below the last tabulated point: extrapolate the last segment.
        dq_dp = (q[-1] - q[-2]) / (p[-1] - p[-2]) if p[-1] != p[-2] else 0.0
        return k * (q[-1] + dq_dp * (p_wh - p[-1])), k * dq_dp


def bottomhole_pressure(q: float, p_res: float, J: float, model: str, p_b: float) -> float:
    """Flowing bottom-hole pressure for oil rate q (inverse of the IPR)."""
    J_pa = J / _PA_PER_BAR
    if model == "linear":
        return p_res - q / J_pa
    pb = min(p_b, p_res)
    q_b = J_pa * (p_res - pb)
    if q <= q_b:
        return p_res - q / J_pa
    c = 1.0 - (q - q_b) * 1.8 / (J_pa * pb)
    if c < 0.0:
        return -1.0                               # beyond the absolute open flow
    x = (-0.2 + math.sqrt(0.04 + 3.2 * c)) / 1.6
    return x * pb


def tubing_wellhead_pressure(
    p_wf: float, m_dot: float, fluid, T: float, D: float, L: float, dz: float,
    roughness: float, n_seg: int,
) -> float | None:
    """Integrate the tubing from the bottom hole up; None if the pressure
    falls below ~1.2 bar before the wellhead (rate not deliverable)."""
    theta = math.asin(max(-1.0, min(1.0, dz / L)))
    area = math.pi * D * D / 4.0
    ds = L / n_seg

    def gradient(P: float) -> float:
        pvt = fluid.pvt(P, T)
        C_L, rho_l, rho_g, mu_l, mu_g = phase_properties_from_pvt(pvt)
        sigma = surface_tension_from_pvt(pvt, float(fluid.api_gravity))
        rho_ns = C_L * rho_l + (1.0 - C_L) * rho_g
        v_m = m_dot / (rho_ns * area)
        return beggs_brill_gradient(v_m, C_L, rho_l, rho_g, mu_l, mu_g, D, roughness, theta,
                                    sigma).dP_dz_total

    P = p_wf
    for _ in range(n_seg):
        g1 = gradient(P)
        P_mid = P - 0.5 * ds * g1
        if P_mid < _P_MIN_PA:
            return None
        P = P - ds * gradient(P_mid)
        if P < _P_MIN_PA:
            return None
    return P


@lru_cache(maxsize=256)
def _curve_cached(api, gg, gor, wor, corr, p_res, T, J, model, D, L, dz, rough, n_seg) -> WellCurve:
    from angelica.properties.black_oil import BlackOilFluid

    fluid = BlackOilFluid(api_gravity=api, gas_gravity=gg, gor_sc_m3_per_m3=gor,
                          wor_sc_m3_per_m3=wor, correlation=corr)
    p_b = fluid.pvt(p_res, T).bubble_point_pa if model == "vogel" else p_res
    m_per = mass_per_oil_sc_m3(fluid)
    J_pa = J / _PA_PER_BAR
    pb = min(p_b, p_res)
    q_aof = (J_pa * p_res if model == "linear"
             else J_pa * (p_res - pb) + J_pa * pb / 1.8)
    rates, pressures = [], []
    for i in range(_N_RATES):
        q = q_aof * (i / (_N_RATES - 1)) ** 1.5   # denser near q = 0, where the peak lies
        p_wf = bottomhole_pressure(q, p_res, J, model, p_b)
        if p_wf <= _P_MIN_PA:
            break
        p_wh = tubing_wellhead_pressure(p_wf, q * m_per / _S_PER_DAY, fluid, T, D, L, dz, rough, n_seg)
        if p_wh is None:
            break
        rates.append(q)
        pressures.append(p_wh)
    if len(rates) < 2:
        rates, pressures = [0.0, 1e-9], [-1.0, -1.0]   # never flows
    peak = max(range(len(pressures)), key=pressures.__getitem__)
    return WellCurve(tuple(rates), tuple(pressures), m_per, peak)


def well_curve(comp, fluid, p_res: float, T: float) -> WellCurve:
    """Deliverability curve of a ``Well`` component for a black-oil fluid."""
    return _curve_cached(
        float(fluid.api_gravity), float(fluid.gas_gravity), float(fluid.gor_sc_m3_per_m3),
        float(fluid.wor_sc_m3_per_m3), str(getattr(fluid, "correlation", "STANDING")),
        round(float(p_res), 0), round(float(T), 3),
        float(comp.productivity_index_m3_per_day_per_bar), comp.ipr_model,
        float(comp.diameter_m), float(comp.tubing_length_m), float(comp.tubing_height_change_m),
        float(comp.tubing_roughness_m), int(comp.tubing_segments),
    )


def well_link_curve(link_state, fluid) -> WellCurve:
    T = link_state.start_node.temperature_c
    if T is None:
        T = getattr(fluid, "reference_temperature_c", 60.0)
    return well_curve(link_state.component, fluid, float(link_state.start_node.pressure_pa), T)
