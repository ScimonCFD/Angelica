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

Tubing temperature
------------------
The fluid enters the bottom hole at the enthalpy it has in the reservoir
(isenthalpic inflow: h(P_wf, T_bh) = h(P_res, T_res)) and climbs the tubing
losing heat to the formation, whose undisturbed temperature falls linearly
from the reservoir temperature (geothermal gradient).  Along the tubing

    dh/ds = −g·sinθ − q'/ṁ,        q' = C·(T − T_e(s))            (first law)
    dT/ds = [dh/ds − (∂h/∂P)_T·dP/ds] / cp

so the fluid temperature also carries the Joule-Thomson effect and the
potential energy of the lift.  The conductance per unit length C is Ramey's
(1962, JPT 14(4), 427): wellbore resistance 1/(2π·r·U) in series with the
transient conduction in the formation f(t)/(2π·k_e),

    C = 2π·r·U·k_e / (k_e + r·U·f(t)),      → 2π·k_e / f(t) for U → ∞

(Ramey's relaxation distance is A = ṁ·cp/C).  f(t) and the expression of A
were read in the OpenGeoSys benchmark documentation ("Wellbore Heat
Transport - EUBHE", which reproduces Ramey's solution — secondary source):

    t_D = α·t / r_w²
    f = 1.1281·√t_D·(1 − 0.3·√t_D)                    t_D ≤ 1.5
    f = (0.4063 + 0.5·ln t_D)·(1 + 0.6/t_D)           t_D > 1.5

For an incompressible fluid in a static column (dh = cp·dT + dP/ρ,
dP = −ρ·g·sinθ·ds) the equations reduce to Ramey's cp·dT/ds = −q'/ṁ, whose
solution for a linear T_e is T(s) = T_e(s) + g_G·A·(1 − e^(−s/A)).
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache

from .beggs_brill import beggs_brill_gradient, phase_properties_from_pvt, surface_tension_from_pvt
from .well_inflow import IPR, ipr_flowing_pressure, ipr_of, ipr_open_flow, mass_per_oil_sc_m3, mass_per_rate_unit

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
    """Flowing bottom-hole pressure for oil rate q, linear or Vogel IPR
    (inverse of the IPR; −1 beyond the absolute open flow).  For any IPR use
    :func:`angelica.closures.well_inflow.ipr_flowing_pressure`."""
    return ipr_flowing_pressure(IPR(model, productivity_index=J), q, p_res, p_b)


_GRAVITY = 9.80665


def ramey_time_function(t_d: float) -> float:
    """Transient heat-conduction function f(t_D) of the formation
    (t_D = α·t/r_w²); see the module docstring for the source."""
    t_d = max(t_d, 1e-12)
    if t_d <= 1.5:
        s = math.sqrt(t_d)
        return 1.1281 * s * (1.0 - 0.3 * s)
    return (0.4063 + 0.5 * math.log(t_d)) * (1.0 + 0.6 / t_d)


def formation_conductance_w_per_mk(
    tubing_diameter_m: float, wellbore_diameter_m: float, u_wellbore_w_per_m2k: float,
    formation_conductivity_w_per_mk: float, formation_diffusivity_m2_per_s: float,
    production_time_s: float,
) -> float:
    """Heat conductance per unit tubing length (W/(m·K)) between the flowing
    fluid and the undisturbed formation: wellbore resistance (``U`` referred
    to the tubing inner radius; inf or nan = none) in series with Ramey's
    transient formation conduction."""
    k_e = formation_conductivity_w_per_mk
    if k_e <= 0.0:
        return 0.0
    r_w = 0.5 * wellbore_diameter_m
    f = ramey_time_function(formation_diffusivity_m2_per_s * production_time_s / (r_w * r_w))
    r = 0.5 * tubing_diameter_m
    U = u_wellbore_w_per_m2k
    if U is None or math.isnan(U) or math.isinf(U):
        return 2.0 * math.pi * k_e / f
    if U <= 0.0:
        return 0.0
    return 2.0 * math.pi * r * U * k_e / (k_e + r * U * f)


@dataclass(frozen=True)
class TubingThermal:
    """Thermal data of a tubing string for :func:`integrate_tubing`.

    ``conductance_w_per_mk`` : C, fluid-to-formation, per unit length
    ``formation_gradient_k_per_m`` : fall of the formation temperature per
        metre of *tubing length* going up (geothermal gradient × dz/L)
    ``formation_bottom_temperature_c`` : formation temperature at the bottom
    """

    conductance_w_per_mk: float
    formation_gradient_k_per_m: float
    formation_bottom_temperature_c: float

    def formation_temperature(self, s: float) -> float:
        return self.formation_bottom_temperature_c - self.formation_gradient_k_per_m * s


@dataclass(frozen=True)
class TubingResult:
    """Outcome of a tubing integration from the bottom hole to the wellhead."""

    wellhead_pressure_pa: float
    wellhead_temperature_c: float
    bottomhole_temperature_c: float
    formation_heat_w: float            # heat given to the formation (> 0: fluid cools)
    wellhead_enthalpy_j_per_kg: float | None
    profile: tuple[tuple[float, float, float, float], ...] = ()   # (s, P, T, T_formation)


def integrate_tubing(
    p_wf: float, t_bh: float, m_dot: float, D: float, L: float, dz: float, n_seg: int,
    gradient, thermal: TubingThermal | None = None, dh_dp=None, cp=None,
    enthalpy=None, temperature_from_h=None, sync_every: int = 0, keep_profile: bool = False,
) -> TubingResult | None:
    """March pressure (midpoint rule) and temperature (exponential
    integrator) up the tubing.

    ``gradient(P, T)`` is the total pressure gradient (Pa/m, positive when the
    pressure falls going up).  With ``thermal`` None the tubing is isothermal.
    Otherwise ``dh_dp(P, T)`` and ``cp(P, T)`` give the temperature change of
    the first law (module docstring).  The enthalpy is tracked exactly,
    h_end = h_start − g·dz − Q/ṁ; when ``enthalpy`` and ``temperature_from_h``
    are given the temperature is re-synchronised with it every ``sync_every``
    segments, which removes the drift of the cp-based estimate (at very low
    rates, where the fluid follows the formation temperature, the heat is
    taken from the enthalpy balance instead).

    Returns None if the pressure falls below ~1.2 bar before the wellhead
    (rate not deliverable).
    """
    ds = L / n_seg
    sin_t = max(-1.0, min(1.0, dz / L))
    P, T = p_wf, t_bh
    q_form = 0.0
    h = enthalpy(P, T) if (thermal is not None and enthalpy is not None) else None
    profile = [(0.0, P, T, thermal.formation_temperature(0.0) if thermal else T)] if keep_profile else []
    m_abs = max(abs(m_dot), 1e-12)
    stiff = False
    for i in range(n_seg):
        g1 = gradient(P, T)
        P_mid = P - 0.5 * ds * g1
        if P_mid < _P_MIN_PA:
            return None
        T_mid = T
        g2 = gradient(P_mid, T_mid)
        if thermal is not None:
            # Over the segment, with θ = T − T_e(s) and T_e falling at γ per metre:
            #   dθ/ds = −θ/A + S + γ,   A = ṁ·cp/C,   S = (−g·sinθ + (∂h/∂P)_T·dP/ds|up) / cp
            # integrated exactly (stable for any ds/A, e.g. at very low rates).
            C = thermal.conductance_w_per_mk
            gamma = thermal.formation_gradient_k_per_m
            theta0 = T - thermal.formation_temperature(i * ds)
            T_new, q_seg = T, 0.0
            for _ in range(2):                          # second pass at the midpoint temperature
                cp_m = cp(P_mid, T_mid)
                S = (-_GRAVITY * sin_t + dh_dp(P_mid, T_mid) * g2) / cp_m
                if C > 0.0:
                    A = m_abs * cp_m / C
                    theta_inf = A * (S + gamma)
                    decay = math.exp(-ds / A)
                    stiff = ds > 3.0 * A
                    theta1 = theta_inf + (theta0 - theta_inf) * decay
                    q_seg = C * (theta_inf * ds + (theta0 - theta_inf) * A * (1.0 - decay))
                else:
                    theta1 = theta0 + (S + gamma) * ds
                    q_seg = 0.0
                T_new = thermal.formation_temperature((i + 1) * ds) + theta1
                T_mid = 0.5 * (T + T_new)
                g2 = gradient(P_mid, T_mid)
            T = T_new
            q_form += q_seg
        P = P - ds * g2
        if P < _P_MIN_PA:
            return None
        if h is not None and sync_every and ((i + 1) % sync_every == 0 or i == n_seg - 1):
            if stiff:
                # Relaxation distance ≪ segment: the temperature is pinned to
                # the formation's and recovering it from the enthalpy is
                # ill-conditioned (Q/ṁ is a difference of large numbers).
                # Keep T; the heat is what the enthalpy balance implies.
                q_form = m_abs * (h - _GRAVITY * sin_t * (i + 1) * ds - enthalpy(P, T))
            else:
                h_now = h - _GRAVITY * sin_t * (i + 1) * ds - q_form / m_abs
                T = temperature_from_h(h_now, P, T)
        if keep_profile:
            s = (i + 1) * ds
            profile.append((s, P, T, thermal.formation_temperature(s) if thermal else T))
    h_end = None if h is None else h - _GRAVITY * dz - q_form / m_abs
    return TubingResult(P, T, t_bh, q_form, h_end, tuple(profile))


def _gradient_function(fluid, m_dot: float, D: float, L: float, dz: float, roughness: float,
                       correlation: str):
    """Total pressure gradient (Pa/m) of a black-oil fluid in the tubing."""
    from .well_correlations import gray_gradient, hagedorn_brown_gradient

    theta = math.asin(max(-1.0, min(1.0, dz / L)))
    area = math.pi * D * D / 4.0

    def gradient(P: float, T: float) -> float:
        pvt = fluid.pvt(P, T)
        C_L, rho_l, rho_g, mu_l, mu_g = phase_properties_from_pvt(pvt)
        sigma = surface_tension_from_pvt(pvt, float(fluid.api_gravity))
        rho_ns = C_L * rho_l + (1.0 - C_L) * rho_g
        v_m = m_dot / (rho_ns * area)
        if correlation == "hagedorn_brown":
            return hagedorn_brown_gradient(v_m, C_L, rho_l, rho_g, mu_l, mu_g, D, roughness, theta,
                                           sigma, P).dP_dz_total
        if correlation == "gray":
            return gray_gradient(v_m, C_L, rho_l, rho_g, mu_l, mu_g, D, roughness, theta,
                                 sigma, P).dP_dz_total
        return beggs_brill_gradient(v_m, C_L, rho_l, rho_g, mu_l, mu_g, D, roughness, theta,
                                    sigma).dP_dz_total

    return gradient


def tubing_wellhead_pressure(
    p_wf: float, m_dot: float, fluid, T: float, D: float, L: float, dz: float,
    roughness: float, n_seg: int, correlation: str = "beggs_brill",
) -> float | None:
    """Wellhead pressure of an *isothermal* tubing at temperature ``T``; None
    if the pressure falls below ~1.2 bar before the wellhead (rate not
    deliverable).  ``correlation``: "beggs_brill", "hagedorn_brown" or "gray".
    The ``Well`` component uses :func:`well_tubing_state`, which also solves
    the temperature."""
    res = integrate_tubing(p_wf, T, m_dot, D, L, dz, n_seg,
                           _gradient_function(fluid, m_dot, D, L, dz, roughness, correlation))
    return None if res is None else res.wellhead_pressure_pa


def _thermal_key(comp) -> tuple:
    """Hashable thermal description of a ``Well``."""
    return ("ramey", float(comp.geothermal_gradient_k_per_m), float(comp.wellbore_diameter_m),
            float(comp.tubing_u_w_per_m2k), float(comp.formation_conductivity_w_per_mk),
            float(comp.formation_diffusivity_m2_per_s), float(comp.production_time_days))


def _black_oil_march(fluid, p_res: float, t_res: float, q_oil: float, ipr: IPR, p_b: float,
                     m_per: float, D: float, L: float, dz: float, rough: float, n_seg: int,
                     tubing_corr: str, thermal_key: tuple, exact: bool,
                     keep_profile: bool = False) -> TubingResult | None:
    """Bottom-hole state from the IPR, then the tubing march, for one oil rate."""
    from angelica.properties.black_oil import _dh_dp, black_oil_enthalpy

    q_ref = q_oil * m_per / mass_per_rate_unit(fluid, ipr)     # oil, or gas for back-pressure
    p_wf = ipr_flowing_pressure(ipr, q_ref, p_res, p_b)
    if p_wf <= _P_MIN_PA:
        return None
    m_dot = q_oil * m_per / _S_PER_DAY
    gradient = _gradient_function(fluid, m_dot, D, L, dz, rough, tubing_corr)
    if thermal_key[0] == "isothermal":
        return integrate_tubing(p_wf, t_res, m_dot, D, L, dz, n_seg, gradient, keep_profile=keep_profile)
    _, g_geo, d_wb, u_wb, k_e, alpha, days = thermal_key
    params = fluid._pvt_params()

    def dh_dp(P: float, T: float) -> float:
        return _dh_dp(params, P, T)

    def cp(P: float, T: float) -> float:
        return fluid.pvt(P, T).mixture_specific_heat_j_per_kg_k

    def enthalpy(P: float, T: float) -> float:
        return black_oil_enthalpy(params, P, T)

    def temperature_from_h(h: float, P: float, T0: float) -> float:
        # two secant steps from the cp-based estimate (already within ~0.1 K)
        h0 = enthalpy(P, T0)
        T1 = T0 + (h - h0) / cp(P, T0)
        if abs(T1 - T0) < 1e-6:
            return T1
        h1 = enthalpy(P, T1)
        if h1 == h0:
            return T1
        T2 = T1 + (h - h1) * (T1 - T0) / (h1 - h0)
        if abs(T2 - T1) < 1e-4:
            return T2
        h2 = enthalpy(P, T2)
        return T2 if h2 == h1 else T2 + (h - h2) * (T2 - T1) / (h2 - h1)

    # isenthalpic inflow: h(P_wf, T_bh) = h(P_res, T_res)
    if exact:
        t_bh = temperature_from_h(enthalpy(p_res, t_res), p_wf, t_res)
    else:
        p_m = 0.5 * (p_res + p_wf)
        t_bh = t_res + dh_dp(p_m, t_res) * (p_res - p_wf) / cp(p_m, t_res)
    thermal = TubingThermal(
        formation_conductance_w_per_mk(D, d_wb, u_wb, k_e, alpha, days * _S_PER_DAY),
        g_geo * dz / L, t_res)
    return integrate_tubing(
        p_wf, t_bh, m_dot, D, L, dz, n_seg, gradient, thermal, dh_dp, cp,
        enthalpy if exact else None, temperature_from_h if exact else None,
        sync_every=5 if exact else 0, keep_profile=keep_profile)


@lru_cache(maxsize=256)
def _curve_cached(api, gg, gor, wor, corr, p_res, T, ipr, D, L, dz, rough, n_seg,
                  tubing_corr="beggs_brill", thermal_key=("isothermal",), gas_cp=2_200.0,
                  gas_k=0.035) -> WellCurve:
    from angelica.properties.black_oil import BlackOilFluid

    fluid = BlackOilFluid(api_gravity=api, gas_gravity=gg, gor_sc_m3_per_m3=gor,
                          wor_sc_m3_per_m3=wor, correlation=corr,
                          gas_specific_heat_j_per_kg_k=gas_cp, gas_thermal_conductivity_w_per_mk=gas_k)
    p_b = fluid.pvt(p_res, T).bubble_point_pa if ipr.uses_bubble_point else p_res
    m_per = mass_per_oil_sc_m3(fluid)
    q_aof = ipr_open_flow(ipr, p_res, p_b) * mass_per_rate_unit(fluid, ipr) / m_per   # oil m³/d
    rates, pressures = [], []
    for i in range(_N_RATES):
        q = q_aof * (i / (_N_RATES - 1)) ** 1.5   # denser near q = 0, where the peak lies
        res = _black_oil_march(fluid, p_res, T, q, ipr, p_b, m_per, D, L, dz, rough, n_seg,
                               tubing_corr, thermal_key, exact=False)
        if res is None:
            break
        rates.append(q)
        pressures.append(res.wellhead_pressure_pa)
    if len(rates) < 2:
        rates, pressures = [0.0, 1e-9], [-1.0, -1.0]   # never flows
    peak = max(range(len(pressures)), key=pressures.__getitem__)
    return WellCurve(tuple(rates), tuple(pressures), m_per, peak)


def _well_args(comp, fluid, p_res: float, T: float) -> tuple:
    return (
        float(fluid.api_gravity), float(fluid.gas_gravity), float(fluid.gor_sc_m3_per_m3),
        float(fluid.wor_sc_m3_per_m3), str(getattr(fluid, "correlation", "STANDING")),
        round(float(p_res), 0), round(float(T), 3),
        ipr_of(comp),
        float(comp.diameter_m), float(comp.tubing_length_m), float(comp.tubing_height_change_m),
        float(comp.tubing_roughness_m), int(comp.tubing_segments),
        str(getattr(comp, "tubing_correlation", "beggs_brill")), _thermal_key(comp),
        float(getattr(fluid, "gas_specific_heat_j_per_kg_k", 2_200.0)),
        float(getattr(fluid, "gas_thermal_conductivity_w_per_mk", 0.035)),
    )


def well_curve(comp, fluid, p_res: float, T: float) -> WellCurve:
    """Deliverability curve of a ``Well`` component for a black-oil fluid."""
    return _curve_cached(*_well_args(comp, fluid, p_res, T))


@lru_cache(maxsize=1024)
def _tubing_state_cached(m_dot, keep_profile, api, gg, gor, wor, corr, p_res, T, ipr, D, L, dz,
                         rough, n_seg, tubing_corr, thermal_key, gas_cp, gas_k) -> TubingResult | None:
    from angelica.properties.black_oil import BlackOilFluid

    fluid = BlackOilFluid(api_gravity=api, gas_gravity=gg, gor_sc_m3_per_m3=gor,
                          wor_sc_m3_per_m3=wor, correlation=corr,
                          gas_specific_heat_j_per_kg_k=gas_cp, gas_thermal_conductivity_w_per_mk=gas_k)
    p_b = fluid.pvt(p_res, T).bubble_point_pa if ipr.uses_bubble_point else p_res
    m_per = mass_per_oil_sc_m3(fluid)
    return _black_oil_march(fluid, p_res, T, m_dot * _S_PER_DAY / m_per, ipr, p_b, m_per,
                            D, L, dz, rough, n_seg, tubing_corr, thermal_key, exact=True,
                            keep_profile=keep_profile)


def well_tubing_state(comp, fluid, p_res: float, T: float, m_dot: float,
                      keep_profile: bool = False) -> TubingResult | None:
    """Tubing solution of a ``Well`` at mass rate ``m_dot`` (kg/s): wellhead
    pressure and temperature, bottom-hole temperature, heat given to the
    formation and (``keep_profile``) the profile (s, P, T, T_formation).
    None when the rate is not deliverable."""
    m_key = float(f"{m_dot:.7g}")
    return _tubing_state_cached(m_key, bool(keep_profile), *_well_args(comp, fluid, p_res, T))


def well_link_curve(link_state, fluid) -> WellCurve:
    T = link_state.start_node.temperature_c
    if T is None:
        T = getattr(fluid, "reference_temperature_c", 60.0)
    return well_curve(link_state.component, fluid, float(link_state.start_node.pressure_pa), T)
