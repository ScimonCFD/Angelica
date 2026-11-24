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

from .beggs_brill import beggs_brill_gradient
from .well_inflow import IPR, ipr_flowing_pressure, ipr_of, ipr_open_flow

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

    @property
    def static_wellhead_pressure_pa(self) -> float:
        """Wellhead pressure of the shut-in (q = 0) column."""
        return self.wellhead_pressures_pa[0]

    @property
    def minimum_stable_rate_kg_per_s(self) -> float:
        """Rate at the maximum of the curve: the smallest rate on the stable branch."""
        return self.oil_rates_sc_m3_per_day[self.peak_index] * self.mass_per_oil_sc_m3 / _S_PER_DAY

    # Hysteresis of a bistable well.  Between the static wellhead pressure and
    # the curve's maximum a well can be dead or producing on the stable
    # branch: a producing well dies when P_wh reaches the maximum; a dead well
    # starts again only when P_wh falls below the static pressure, where the
    # shut-in column itself would push fluid out.
    def dies_at(self, p_wh: float) -> bool:
        return p_wh >= self.max_wellhead_pressure_pa

    def starts_at(self, p_wh: float) -> bool:
        return p_wh < self.static_wellhead_pressure_pa

    def mass_flow(self, p_wh: float, dead: bool = False) -> tuple[float, float]:
        """(ṁ, dṁ/dP_wh) on the stable branch at wellhead pressure ``p_wh``;
        (0, 0) above the maximum, or when the well is ``dead`` (see
        :meth:`dies_at` / :meth:`starts_at`)."""
        q, p = self.oil_rates_sc_m3_per_day, self.wellhead_pressures_pa
        k = self.mass_per_oil_sc_m3 / _S_PER_DAY
        i0 = self.peak_index
        if p_wh >= p[i0] or dead:
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
    wellhead_enthalpy_j_per_kg: float | None   # of the wellhead stream (with lift gas)
    profile: tuple[tuple[float, float, float, float], ...] = ()   # (s, P, T, T_formation)
    # gas lift (zero / nan without it)
    lift_mass_flow_kg_per_s: float = 0.0
    injection_pressure_pa: float = float("nan")
    injection_temperature_c: float = float("nan")      # lift gas arriving at the valve
    outlet_fluid: object = None                       # adapter of the wellhead stream


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


def _gradient_function(wf, m_dot: float, D: float, L: float, dz: float, roughness: float,
                       correlation: str):
    """Total pressure gradient (Pa/m) in the tubing; ``wf`` is a well-fluid
    adapter (:mod:`angelica.closures.well_fluid`) or a black-oil fluid."""
    from .well_correlations import gray_gradient, hagedorn_brown_gradient
    from .well_fluid import BlackOilWellFluid

    if hasattr(wf, "pvt"):
        wf = BlackOilWellFluid(wf)
    theta = math.asin(max(-1.0, min(1.0, dz / L)))
    area = math.pi * D * D / 4.0

    def gradient(P: float, T: float) -> float:
        C_L, rho_l, rho_g, mu_l, mu_g, sigma = wf.phase_properties(P, T)
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
    ``fluid``: a black-oil fluid or a well-fluid adapter.  The ``Well``
    component uses :func:`well_tubing_state`, which also solves the
    temperature."""
    res = integrate_tubing(p_wf, T, m_dot, D, L, dz, n_seg,
                           _gradient_function(fluid, m_dot, D, L, dz, roughness, correlation))
    return None if res is None else res.wellhead_pressure_pa


def _thermal_key(comp) -> tuple:
    """Hashable thermal description of a ``Well``."""
    return ("ramey", float(comp.geothermal_gradient_k_per_m), float(comp.wellbore_diameter_m),
            float(comp.tubing_u_w_per_m2k), float(comp.formation_conductivity_w_per_mk),
            float(comp.formation_diffusivity_m2_per_s), float(comp.production_time_days))


def _march(wf, p_res: float, t_res: float, q: float, ipr: IPR, p_b: float,
           m_per: float, D: float, L: float, dz: float, rough: float, n_seg: int,
           tubing_corr: str, thermal_key: tuple, exact: bool,
           keep_profile: bool = False, lift: tuple | None = None) -> TubingResult | None:
    """Bottom-hole state from the IPR, then the tubing march, for one rate
    ``q`` of the curve's base unit (``m_per`` kg per standard m³).

    ``lift`` = (lift-gas adapter, lift mass rate kg/s, injection depth m
    below the wellhead along the tubing) for a gas-lifted well."""
    q_ref = q * m_per / wf.mass_per_rate(ipr)     # rate of the IPR's reference phase
    p_wf = ipr_flowing_pressure(ipr, q_ref, p_res, p_b)
    if p_wf <= _P_MIN_PA:
        return None
    m_dot = q * m_per / _S_PER_DAY
    thermal = None
    if thermal_key[0] == "isothermal":
        t_bh = t_res
    else:
        _, g_geo, d_wb, u_wb, k_e, alpha, days = thermal_key
        # isenthalpic inflow: h(P_wf, T_bh) = h(P_res, T_res)
        if exact:
            t_bh = wf.temperature_from_h(wf.enthalpy(p_res, t_res), p_wf, t_res)
        else:
            p_m = 0.5 * (p_res + p_wf)
            t_bh = t_res + wf.dh_dp(p_m, t_res) * (p_res - p_wf) / wf.cp(p_m, t_res)
        thermal = TubingThermal(
            formation_conductance_w_per_mk(D, d_wb, u_wb, k_e, alpha, days * _S_PER_DAY),
            g_geo * dz / L, t_res)

    def section(fl, p0, t0, m, length, rise, n, therm, keep):
        gradient = _gradient_function(fl, m, D, length, rise, rough, tubing_corr)
        if therm is None:
            return integrate_tubing(p0, t0, m, D, length, rise, n, gradient, keep_profile=keep)
        return integrate_tubing(
            p0, t0, m, D, length, rise, n, gradient, therm, fl.dh_dp, fl.cp,
            fl.enthalpy if exact else None, fl.temperature_from_h if exact else None,
            sync_every=5 if exact else 0, keep_profile=keep)

    if lift is None:
        return section(wf, p_wf, t_bh, m_dot, L, dz, n_seg, thermal, keep_profile)

    from .well_fluid import GasLiftMixture

    gas, m_gl, depth = lift
    s_inj = L - depth                                    # injection point, from the bottom
    n1 = int(round(n_seg * s_inj / L))
    profile: list = []
    q_form = 0.0
    if n1 >= 1:
        low = section(wf, p_wf, t_bh, m_dot, s_inj, dz * s_inj / L, n1, thermal, keep_profile)
        if low is None:
            return None
        p1, t1, q_form = low.wellhead_pressure_pa, low.wellhead_temperature_c, low.formation_heat_w
        profile = list(low.profile)
    else:
        p1, t1, s_inj = p_wf, t_bh, 0.0
    # the lift gas arrives at the formation temperature of the injection depth
    t_gas = thermal.formation_temperature(s_inj) if thermal is not None else t1
    m_tot = m_dot + m_gl
    mix = GasLiftMixture(wf, gas, m_gl / m_tot)
    h_mix = ((m_dot * wf.enthalpy(p1, t1) if m_dot > 0.0 else 0.0) + m_gl * gas.enthalpy(p1, t_gas)) / m_tot
    t2 = mix.temperature_from_h(h_mix, p1, t1 if m_dot > 0.0 else t_gas)
    therm2 = None if thermal is None else TubingThermal(
        thermal.conductance_w_per_mk, thermal.formation_gradient_k_per_m, thermal.formation_temperature(s_inj))
    up = section(mix, p1, t2, m_tot, L - s_inj, dz * (L - s_inj) / L, max(n_seg - n1, 1), therm2, keep_profile)
    if up is None:
        return None
    if keep_profile:
        profile += [(s_inj + s_, P_, T_, Te_) for (s_, P_, T_, Te_) in up.profile[(1 if profile else 0):]]
    return TubingResult(up.wellhead_pressure_pa, up.wellhead_temperature_c, t_bh, q_form + up.formation_heat_w,
                        up.wellhead_enthalpy_j_per_kg, tuple(profile), m_gl, p1, t_gas, mix)


def _bubble_point_for(wf, ipr: IPR, p_res: float, T: float) -> float:
    return wf.bubble_point_pa(p_res, T) if ipr.uses_bubble_point else p_res


@lru_cache(maxsize=256)
def _curve_cached(wf, p_res, T, ipr, D, L, dz, rough, n_seg, tubing_corr="beggs_brill",
                  thermal_key=("isothermal",), lift=None) -> WellCurve:
    p_b = _bubble_point_for(wf, ipr, p_res, T)
    m_per = wf.mass_per_base_unit_for(ipr)
    q_aof = ipr_open_flow(ipr, p_res, p_b) * wf.mass_per_rate(ipr) / m_per   # base units, m³/d
    rates, pressures = [], []
    for i in range(_N_RATES):
        q = q_aof * (i / (_N_RATES - 1)) ** 1.5   # denser near q = 0, where the peak lies
        res = _march(wf, p_res, T, q, ipr, p_b, m_per, D, L, dz, rough, n_seg,
                     tubing_corr, thermal_key, exact=False, lift=lift)
        if res is None:
            break
        rates.append(q)
        pressures.append(res.wellhead_pressure_pa)
    if len(rates) < 2:
        rates, pressures = [0.0, 1e-9], [-1.0, -1.0]   # never flows
    peak = max(range(len(pressures)), key=pressures.__getitem__)
    return WellCurve(tuple(rates), tuple(pressures), m_per, peak)


def _well_args(comp, fluid, p_res: float, T: float, link_state=None) -> tuple:
    from .well_fluid import well_fluid_for

    p_res, T = round(float(p_res), 0), round(float(T), 3)
    wf = well_fluid_for(fluid, link_state, comp, p_res, T)
    if wf is None:
        raise ValueError("Well links need a black-oil or compositional fluid.")
    return (
        wf, p_res, T, ipr_of(comp),
        float(comp.diameter_m), float(comp.tubing_length_m), float(comp.tubing_height_change_m),
        float(comp.tubing_roughness_m), int(comp.tubing_segments),
        str(getattr(comp, "tubing_correlation", "beggs_brill")), _thermal_key(comp),
        _lift_args(comp, wf),
    )


def _lift_args(comp, wf) -> tuple | None:
    """(lift-gas adapter, lift mass rate kg/s, injection depth m) or None."""
    if not getattr(comp, "has_gas_lift", False):
        return None
    from .well_fluid import lift_gas_fluid, lift_gas_mass_per_sc_m3

    gas = lift_gas_fluid(wf, comp)
    m_gl = float(comp.gas_lift_rate_sm3_per_day) * lift_gas_mass_per_sc_m3(gas) / _S_PER_DAY
    return (gas, round(m_gl, 12), float(comp.gas_lift_depth_m))


def gas_lift_mass_flow(comp, fluid, p_res: float, T: float, link_state=None) -> float:
    """Lift-gas mass rate (kg/s) of a ``Well``; 0 without gas lift."""
    lift = _well_args(comp, fluid, p_res, T, link_state)[-1]
    return 0.0 if lift is None else lift[1]


def well_curve(comp, fluid, p_res: float, T: float, link_state=None) -> WellCurve:
    """Deliverability curve of a ``Well`` component.  ``fluid``: black-oil,
    or compositional (the composition is the link's, ``link_state.zs``, or
    the fluid's default).  For a compositional well the curve's rates are
    of the IPR's reference phase (gas for the back-pressure law)."""
    return _curve_cached(*_well_args(comp, fluid, p_res, T, link_state))


@lru_cache(maxsize=1024)
def _tubing_state_cached(m_dot, keep_profile, wf, p_res, T, ipr, D, L, dz, rough, n_seg,
                         tubing_corr, thermal_key, lift=None) -> TubingResult | None:
    p_b = _bubble_point_for(wf, ipr, p_res, T)
    m_per = wf.mass_per_base_unit_for(ipr)
    return _march(wf, p_res, T, m_dot * _S_PER_DAY / m_per, ipr, p_b, m_per,
                  D, L, dz, rough, n_seg, tubing_corr, thermal_key, exact=True,
                  keep_profile=keep_profile, lift=lift)


def well_tubing_state(comp, fluid, p_res: float, T: float, m_dot: float,
                      keep_profile: bool = False, link_state=None) -> TubingResult | None:
    """Tubing solution of a ``Well`` at mass rate ``m_dot`` (kg/s): wellhead
    pressure and temperature, bottom-hole temperature, heat given to the
    formation and (``keep_profile``) the profile (s, P, T, T_formation).
    None when the rate is not deliverable."""
    m_key = float(f"{m_dot:.7g}")
    return _tubing_state_cached(m_key, bool(keep_profile), *_well_args(comp, fluid, p_res, T, link_state))


def well_link_curve(link_state, fluid) -> WellCurve:
    T = link_state.start_node.temperature_c
    if T is None:
        T = getattr(fluid, "reference_temperature_c", 60.0)
    return well_curve(link_state.component, fluid, float(link_state.start_node.pressure_pa), T, link_state)
