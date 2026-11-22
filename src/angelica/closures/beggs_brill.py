"""
Beggs-Brill (1973) two-phase pipe-flow pressure-gradient correlation.

Reference
---------
Beggs, H.D., Brill, J.P. (1973). A Study of Two-Phase Flow in Inclined
Pipes. Journal of Petroleum Technology, May, 607-617.

Brill, J.P., Mukherjee, H. (1999). Multiphase Flow in Wells. SPE
Monograph 17, §4.2.3 (flow-pattern map with the transition region of the
revised correlation, and the equations used here).

Inclination correction (original Beggs-Brill)
---------------------------------------------
    H_L(θ) = H_L(0)·ψ,   ψ = 1 + C·[sin(1.8θ) − ⅓·sin³(1.8θ)]
    C = (1 − λ_L)·ln(e·λ_L^f·N_LV^g·N_Fr^h),   C ≥ 0
    N_LV = v_sL·(ρ_L / (g·σ_L))^¼   (liquid velocity number)

    uphill   segregated    e = 0.011  f = −3.768   g =  3.539   h = −1.614
             intermittent  e = 2.96   f =  0.305   g = −0.4473  h =  0.0978
             distributed   C = 0
    downhill all patterns  e = 4.70   f = −0.3692  g =  0.1244  h = −0.5056

θ is the angle from horizontal, negative downhill (ψ < 1: downhill flow
holds less liquid).  In the transition region the segregated and
intermittent holdups are each corrected, then weighted.

Payne correction (optional)
---------------------------
Payne, W.R., Palmer, C.M., Brill, J.P. (1979). Evaluation of Inclined
Pipe, Two-Phase Liquid Holdup and Pressure-Drop Correlations Using
Experimental Data. JPT, Sept., 1198-1208: Beggs-Brill over-predicts
holdup, so H_L(θ) is multiplied by 0.924 uphill and 0.685 downhill.
Enabled with ``BeggsBrillCorrelation(payne_correction=True)``.

Physics
-------
The no-slip model (homogeneous flow) that Angelica uses by default
assumes that gas and liquid travel at the same velocity.  In reality,
for upward and horizontal flow the gas travels faster than the liquid
(slip), so the in-situ liquid fraction (holdup H_L) is greater than
the input liquid content C_L.  Beggs-Brill quantifies that holdup and
uses it to correct:
  • the gravity term  (uses in-situ density ρ_s = H_L·ρ_l + (1-H_L)·ρ_g)
  • the friction term (two-phase friction factor f_tp ≥ f_n)

Integration
-----------
`BeggsBrillCorrelation` implements `PressureDropCorrelation`.  Because
the interface does not pass the fluid model to `calculate_velocity`, the
correlation must be told which fluid it will see via `set_fluid()` before
the first solve.  `SteadyBlackOilSolver` calls this automatically when
the turbulent correlation is an instance of `BeggsBrillCorrelation`.
"""
from __future__ import annotations

import math
import warnings
from dataclasses import dataclass

from .gravity import GRAVITY_M_PER_S2
from .pressure_drop import PressureDropCorrelation

# ── Constants ──────────────────────────────────────────────────────────────────

_G = GRAVITY_M_PER_S2       # m/s²
_MIN_VELOCITY = 1e-6        # m/s — below this treat flow as zero

# ── Standalone friction-factor helper (Churchill 1977) ────────────────────────

def _darcy_friction_factor(re: float, rel_roughness: float) -> float:
    """Darcy friction factor for any Re, explicit Churchill (1977) formula."""
    re = max(re, 1.0)
    if re < 2300.0:
        return 64.0 / re
    A = (2.457 * math.log(1.0 / ((7.0 / re) ** 0.9 + 0.27 * rel_roughness))) ** 16
    B = (37530.0 / re) ** 16
    return 8.0 * ((8.0 / re) ** 12 + (A + B) ** (-1.5)) ** (1.0 / 12.0)


# ── Flow regime ───────────────────────────────────────────────────────────────

def _regime_limits(C_L: float) -> tuple[float, float, float, float]:
    """Beggs-Brill flow-pattern transition boundaries (Froude number)."""
    L1 = 316.0   * max(C_L, 1e-9) ** 0.302
    L2 = 9.252e-4 * max(C_L, 1e-9) ** (-2.4684)
    L3 = 0.10    * max(C_L, 1e-9) ** (-1.4516)
    L4 = 0.50    * max(C_L, 1e-9) ** (-6.738)
    return L1, L2, L3, L4


def flow_regime(C_L: float, N_Fr: float) -> tuple[str, float]:
    """Identify horizontal flow regime and transition interpolation factor.

    Parameters
    ----------
    C_L : input liquid content (no-slip holdup), 0 ≤ C_L ≤ 1
    N_Fr: mixture Froude number = v_m² / (g·D)

    Returns
    -------
    (regime, A) where regime ∈ {"segregated","transition","intermittent","distributed"}
    and A ∈ [0,1] is used only in the transition regime as a linear weight
    between segregated (A=1) and intermittent (A=0).
    """
    if C_L <= 0.0:
        return "distributed", 0.0
    if C_L >= 1.0:
        return "intermittent", 0.0

    L1, L2, L3, L4 = _regime_limits(C_L)

    if (C_L < 0.01 and N_Fr < L1) or (C_L >= 0.01 and N_Fr < L2):
        return "segregated", 0.0
    if C_L >= 0.01 and L2 <= N_Fr < L3:
        span = L3 - L2
        A = (L3 - N_Fr) / span if span > 0.0 else 0.5
        return "transition", float(A)
    if (C_L < 0.4 and L3 <= N_Fr <= L1) or (C_L >= 0.4 and L3 <= N_Fr <= L4):
        return "intermittent", 0.0
    return "distributed", 0.0


# ── Liquid holdup ─────────────────────────────────────────────────────────────

_HOLDUP_COEFFS = {
    "segregated":   (0.98,  0.4846, 0.0868),
    "intermittent": (0.845, 0.5351, 0.0173),
    "distributed":  (1.065, 0.5824, 0.0609),
}
_INCLINATION_UPHILL = {
    "segregated":   (0.011, -3.768, 3.539, -1.614),
    "intermittent": (2.96, 0.305, -0.4473, 0.0978),
}
_INCLINATION_DOWNHILL = (4.70, -0.3692, 0.1244, -0.5056)
_PAYNE_UPHILL = 0.924
_PAYNE_DOWNHILL = 0.685


def _holdup_horizontal_pattern(pattern: str, C_L: float, N_Fr: float) -> float:
    a, b, c = _HOLDUP_COEFFS[pattern]
    HL0 = a * C_L ** b / max(N_Fr ** c, 1e-30)
    return max(min(HL0, 1.0), C_L)          # H_L(0) ≥ λ_L


def holdup_horizontal(C_L: float, N_Fr: float) -> float:
    """Horizontal liquid holdup H_L(0) (Beggs-Brill 1973), H_L(0) ∈ [C_L, 1]."""
    if C_L <= 0.0:
        return 0.0
    if C_L >= 1.0:
        return 1.0
    regime, A = flow_regime(C_L, N_Fr)
    if regime == "transition":
        return (A * _holdup_horizontal_pattern("segregated", C_L, N_Fr)
                + (1.0 - A) * _holdup_horizontal_pattern("intermittent", C_L, N_Fr))
    return _holdup_horizontal_pattern(regime, C_L, N_Fr)


def inclination_factor(pattern: str, C_L: float, N_Fr: float, N_LV: float, theta_rad: float) -> float:
    """Beggs-Brill inclination factor ψ for one flow pattern (θ > 0 uphill)."""
    if abs(theta_rad) < 1e-9 or C_L <= 0.0 or C_L >= 1.0:
        return 1.0
    if theta_rad > 0.0:
        if pattern not in _INCLINATION_UPHILL:
            return 1.0                          # distributed uphill: C = 0
        e, f, g, h = _INCLINATION_UPHILL[pattern]
    else:
        e, f, g, h = _INCLINATION_DOWNHILL
    arg = e * C_L ** f * max(N_LV, 1e-12) ** g * max(N_Fr, 1e-12) ** h
    C = (1.0 - C_L) * math.log(arg) if arg > 0.0 else 0.0
    C = max(C, 0.0)
    s18 = math.sin(1.8 * theta_rad)
    return 1.0 + C * (s18 - s18 ** 3 / 3.0)


def liquid_holdup(
    C_L: float, N_Fr: float, N_LV: float, theta_rad: float, payne_correction: bool = False,
) -> float:
    """In-situ liquid holdup H_L(θ): Beggs-Brill with inclination correction."""
    if C_L <= 0.0:
        return 0.0
    if C_L >= 1.0:
        return 1.0
    regime, A = flow_regime(C_L, N_Fr)

    def _pattern(p: str) -> float:
        return (_holdup_horizontal_pattern(p, C_L, N_Fr)
                * inclination_factor(p, C_L, N_Fr, N_LV, theta_rad))

    if regime == "transition":
        HL = A * _pattern("segregated") + (1.0 - A) * _pattern("intermittent")
    else:
        HL = _pattern(regime)
    if payne_correction and abs(theta_rad) > 1e-9:
        HL *= _PAYNE_UPHILL if theta_rad > 0.0 else _PAYNE_DOWNHILL
    return max(min(HL, 1.0), 1e-9)


def liquid_velocity_number(v_sl: float, rho_l: float, sigma_l: float) -> float:
    """N_LV = v_sL·(ρ_L/(g·σ_L))^¼ (dimensionless, SI units)."""
    return abs(v_sl) * (rho_l / (_G * max(sigma_l, 1e-6))) ** 0.25


# ── Two-phase friction factor ─────────────────────────────────────────────────

def two_phase_friction_factor(f_n: float, C_L: float, H_L: float) -> float:
    """Beggs-Brill two-phase friction factor multiplier.

    Parameters
    ----------
    f_n : no-slip (single-phase) Darcy friction factor at Re_ns
    C_L : input liquid content
    H_L : in-situ liquid holdup (after inclination correction)

    Returns f_tp = f_n * exp(s), where s is the Beggs-Brill friction
    correction.  If H_L ≈ C_L (no-slip limit) then s → 0 and f_tp → f_n.
    """
    if H_L <= 0.0 or H_L >= 1.0 or C_L <= 0.0:
        return f_n
    y = C_L / max(H_L ** 2, 1e-12)
    if y <= 1.0:
        return f_n
    if y < 1.2:
        s = math.log(max(2.2 * y - 1.2, 1e-30))
    else:
        ln_y = math.log(y)
        denom = -0.0523 + 3.182 * ln_y - 0.8725 * ln_y ** 2 + 0.01853 * ln_y ** 4
        if abs(denom) < 1e-12:
            s = 0.0
        else:
            s = ln_y / denom
    return f_n * math.exp(min(s, 2.0))   # cap at exp(2) ≈ 7.4 to prevent runaway


# ── Pressure gradient ─────────────────────────────────────────────────────────

@dataclass(frozen=True)
class BeggsBrillState:
    """Intermediate Beggs-Brill calculation state (for diagnostics/testing)."""
    regime: str
    holdup: float          # H_L — in-situ liquid holdup
    C_L: float             # no-slip liquid content
    N_Fr: float            # mixture Froude number
    N_LV: float            # liquid velocity number
    rho_ns: float          # no-slip density (kg/m³)
    rho_s: float           # in-situ (slip) density (kg/m³)
    f_n: float             # no-slip friction factor
    f_tp: float            # two-phase friction factor
    dP_dz_friction: float  # Pa/m
    dP_dz_gravity: float   # Pa/m
    dP_dz_total: float     # Pa/m


def beggs_brill_gradient(
    v_m: float,
    C_L: float,
    rho_l: float,
    rho_g: float,
    mu_l: float,
    mu_g: float,
    D: float,
    roughness: float,
    theta_rad: float,
    sigma_l: float = 0.030,
    payne_correction: bool = False,
) -> BeggsBrillState:
    """Compute the Beggs-Brill pressure gradient at a given mixture velocity.

    Parameters
    ----------
    v_m       : mixture velocity (m/s), always positive
    C_L       : no-slip liquid input content (–)
    rho_l     : liquid density (kg/m³)
    rho_g     : gas density (kg/m³)
    mu_l      : liquid dynamic viscosity (Pa·s)
    mu_g      : gas dynamic viscosity (Pa·s)
    D         : pipe inner diameter (m)
    roughness : absolute wall roughness (m)
    theta_rad : pipe inclination angle (rad), positive = uphill
    sigma_l   : gas-liquid surface tension (N/m), for N_LV
    payne_correction : multiply H_L by 0.924 uphill / 0.685 downhill

    Returns
    -------
    BeggsBrillState with all intermediate quantities and Pa/m gradients.
    """
    v_m = max(v_m, _MIN_VELOCITY)

    # ── No-slip mixture properties ─────────────────────────────────────────
    rho_ns = C_L * rho_l + (1.0 - C_L) * rho_g
    mu_ns  = C_L * mu_l  + (1.0 - C_L) * mu_g

    # ── Flow regime and holdup ─────────────────────────────────────────────
    N_Fr = v_m ** 2 / (_G * max(D, 1e-6))
    N_LV = liquid_velocity_number(C_L * v_m, rho_l, sigma_l)
    HL   = liquid_holdup(C_L, N_Fr, N_LV, theta_rad, payne_correction)

    regime, _ = flow_regime(C_L, N_Fr)

    # ── In-situ (slip) density ─────────────────────────────────────────────
    rho_s = HL * rho_l + (1.0 - HL) * rho_g

    # ── No-slip friction factor ────────────────────────────────────────────
    Re_ns = rho_ns * v_m * D / max(mu_ns, 1e-12)
    f_n   = _darcy_friction_factor(Re_ns, roughness / max(D, 1e-6))

    # ── Two-phase friction factor ──────────────────────────────────────────
    f_tp  = two_phase_friction_factor(f_n, C_L, HL)

    # ── Pressure gradients ─────────────────────────────────────────────────
    dP_f = f_tp * rho_ns * v_m ** 2 / (2.0 * D)
    dP_g = rho_s * _G * math.sin(theta_rad)
    dP_total = dP_f + dP_g

    return BeggsBrillState(
        regime=regime, holdup=HL, C_L=C_L, N_Fr=N_Fr, N_LV=N_LV,
        rho_ns=rho_ns, rho_s=rho_s,
        f_n=f_n, f_tp=f_tp,
        dP_dz_friction=dP_f,
        dP_dz_gravity=dP_g,
        dP_dz_total=dP_total,
    )


def phase_properties_from_pvt(pvt) -> tuple[float, float, float, float, float]:
    """(C_L, rho_l, rho_g, mu_l, mu_g) for Beggs-Brill from a black-oil PVT state.

    Liquid = oil + water (volume-weighted); C_L is the no-slip liquid content.
    """
    alpha_l = pvt.holdup_oil + pvt.holdup_water
    alpha_g = pvt.holdup_gas

    if alpha_l + alpha_g < 1e-9:
        return 1.0, pvt.mixture_density_kg_per_m3, 1.2, pvt.mixture_viscosity_pa_s, 1.5e-5

    C_L = alpha_l / (alpha_l + alpha_g)

    # Liquid phase: volume-weighted oil + water
    if alpha_l > 1e-9:
        rho_l = (pvt.holdup_oil * pvt.density_oil_kg_per_m3
                 + pvt.holdup_water * pvt.density_water_kg_per_m3) / alpha_l
        mu_l  = (pvt.holdup_oil * pvt.viscosity_oil_pa_s
                 + pvt.holdup_water * pvt.viscosity_water_pa_s) / alpha_l
    else:
        rho_l = pvt.mixture_density_kg_per_m3
        mu_l  = pvt.mixture_viscosity_pa_s

    rho_g = max(pvt.density_gas_kg_per_m3, 0.1)
    mu_g  = max(pvt.viscosity_gas_pa_s, 1e-6)

    return C_L, rho_l, rho_g, mu_l, mu_g


def surface_tension_from_pvt(pvt, api_gravity: float) -> float:
    """Liquid surface tension (N/m) for a black-oil PVT state: Baker-Swerdloff
    oil and gas-water values, no-slip volume average."""
    from angelica.properties.surface_tension import (
        liquid_surface_tension,
        oil_surface_tension,
        water_surface_tension,
    )
    return liquid_surface_tension(
        oil_surface_tension(api_gravity, pvt.temperature_c, pvt.pressure_pa),
        water_surface_tension(pvt.temperature_c, pvt.pressure_pa),
        pvt.holdup_oil, pvt.holdup_water,
    )


def link_phase_properties(fluid, link_state) -> tuple[float, float, float, float, float, float] | None:
    """``(C_L, rho_l, rho_g, mu_l, mu_g, sigma_l)`` of a link for slip
    correlations, from a black-oil (``.pvt``) or compositional
    (``.beggs_brill_phase_properties``) fluid; None for other fluids."""
    if fluid is not None and hasattr(fluid, "fluid_for_link"):
        fluid = fluid.fluid_for_link(link_state)  # per-inlet composition wrapper
    if fluid is not None and hasattr(fluid, "beggs_brill_phase_properties"):
        return fluid.beggs_brill_phase_properties(link_state)
    if fluid is None or not hasattr(fluid, "pvt"):
        return None
    P_s = getattr(getattr(link_state, "start_node", None), "pressure_pa", None)
    P_e = getattr(getattr(link_state, "end_node", None), "pressure_pa", None)
    if P_s is not None and P_e is not None:
        P_avg = 0.5 * (float(P_s) + float(P_e))
    elif P_s is not None:
        P_avg = float(P_s)
    elif P_e is not None:
        P_avg = float(P_e)
    else:
        P_avg = getattr(fluid, "reference_pressure_pa", 101_325.0)
    T_avg = float(link_state.temperature_c) if link_state.temperature_c is not None else 20.0
    pvt = fluid.pvt(P_avg, T_avg)
    return (*phase_properties_from_pvt(pvt),
            surface_tension_from_pvt(pvt, float(getattr(fluid, "api_gravity", 35.0))))


def multiphase_flow_report(link_state, fluid, correlation=None):
    """Two-phase flow state of a solved pipe (``MultiphaseFlowResult``), or None.

    With a ``BeggsBrillCorrelation`` the holdup, pattern and gradients are
    Beggs-Brill's; otherwise the homogeneous model is reported (H_L = C_L,
    friction = total gradient minus elevation).
    """
    from angelica.core.results import MultiphaseFlowResult

    comp = getattr(link_state, "component", None)
    L = getattr(comp, "length_m", None)
    if not L or not hasattr(comp, "absolute_roughness_m"):
        return None
    props = link_phase_properties(fluid, link_state)
    if props is None:
        return None
    C_L, rho_l, rho_g, mu_l, mu_g, sigma_l = props
    D = comp.diameter_m
    area = math.pi * D * D / 4.0
    m = float(link_state.mass_flow_kg_per_s)
    rho_ns = C_L * rho_l + (1.0 - C_L) * rho_g
    v_m = abs(m) / max(rho_ns * area, 1e-30)
    direction = 1.0 if m >= 0.0 else -1.0
    theta = direction * math.asin(max(-1.0, min(1.0, comp.height_change_m / L)))
    p_up, p_down = float(link_state.start_node.pressure_pa), float(link_state.end_node.pressure_pa)
    if direction < 0.0:
        p_up, p_down = p_down, p_up
    if isinstance(correlation, BeggsBrillCorrelation):
        st = beggs_brill_gradient(v_m, C_L, rho_l, rho_g, mu_l, mu_g, D, comp.absolute_roughness_m,
                                  theta, sigma_l, correlation.payne_correction)
        pattern, H_L = st.regime, st.holdup
        elev, fric = st.dP_dz_gravity, st.dP_dz_friction
    else:
        pattern, H_L = "no-slip", C_L
        elev = rho_ns * _G * math.sin(theta)
        fric = (p_up - p_down) / L - elev
    return MultiphaseFlowResult(
        flow_pattern=pattern,
        no_slip_liquid_fraction=C_L,
        liquid_holdup=H_L,
        superficial_liquid_velocity_m_per_s=C_L * v_m,
        superficial_gas_velocity_m_per_s=(1.0 - C_L) * v_m,
        mixture_velocity_m_per_s=v_m,
        liquid_volume_m3=H_L * area * L,
        elevation_gradient_pa_per_m=elev,
        friction_gradient_pa_per_m=fric,
    )


# ── PressureDropCorrelation implementation ────────────────────────────────────

class BeggsBrillCorrelation(PressureDropCorrelation):
    """Beggs-Brill (1973) two-phase pipe-flow correlation.

    This correlation replaces the single-phase Colebrook-White model for
    pipelines carrying a gas-liquid mixture.  It accounts for slip between
    phases and correctly calculates the in-situ liquid holdup, which affects
    both the gravity and friction pressure gradients.

    Usage
    -----
    Pass as ``turbulent_pipe_correlation`` to ``SteadyBlackOilSolver``.
    The solver automatically calls ``set_fluid()`` before each outer
    iteration to keep the phase-fraction data current::

        from angelica.closures.beggs_brill import BeggsBrillCorrelation
        result = SteadyBlackOilSolver(
            turbulent_pipe_correlation=BeggsBrillCorrelation()
        ).solve(case)

    Limitations
    -----------
    - Requires a ``BlackOilFluid`` (or compatible fluid with ``.pvt(P, T)``).
      Falls back to no-slip single-phase Colebrook-White otherwise.
    - Acceleration pressure drop (small, <1% for most pipelines) is omitted.
    - Friction uses the Churchill (1977) form of the rough-pipe Moody chart
      for the no-slip friction factor, as in the "BBM" variant.

    Parameters
    ----------
    payne_correction : multiply the holdup by 0.924 uphill and 0.685
        downhill (Payne et al. 1979).  Off by default (original correlation).
    """

    def __init__(self, payne_correction: bool = False) -> None:
        self._fluid = None
        self._warned_fallback = False
        self.payne_correction = bool(payne_correction)

    def set_fluid(self, fluid) -> None:
        """Update the fluid model used to extract phase-fraction data."""
        self._fluid = fluid

    def _phase_properties(self, pipe_state) -> tuple[float, float, float, float, float, float]:
        """Return (C_L, rho_l, rho_g, mu_l, mu_g, sigma_l) at midpoint conditions."""
        props = link_phase_properties(self._fluid, pipe_state)
        if props is None:
            if not self._warned_fallback:
                warnings.warn(
                    "BeggsBrillCorrelation: fluid has no black-oil PVT; using an "
                    "all-liquid fallback (C_L=1, rho_l=850 kg/m3).",
                    stacklevel=2,
                )
                self._warned_fallback = True
            return 1.0, 850.0, 1.2, 1e-3, 1.5e-5, 0.030  # all-liquid fallback
        return props

    def _solve_vm(
        self,
        net_driving: float,
        C_L: float,
        rho_l: float,
        rho_g: float,
        mu_l: float,
        mu_g: float,
        D: float,
        L: float,
        roughness: float,
        theta_rad: float,
        v0: float,
        sigma_l: float = 0.030,
    ) -> tuple[float, float, float]:
        """Fixed-point iteration to find v_m such that ΔP_friction = net_driving.

        From Darcy-Weisbach:  ΔP_f = f_tp · ρ_ns · v_m² · L / (2·D)
        Rearranging:          v_m  = √(2·D·ΔP_f / (f_tp · ρ_ns · L))

        Returns (v_m, f_tp, H_L).
        """
        v_m = max(abs(v0), _MIN_VELOCITY)
        f_tp_last = 0.02
        HL_last   = C_L

        for _ in range(60):
            st = beggs_brill_gradient(v_m, C_L, rho_l, rho_g, mu_l, mu_g,
                                      D, roughness, theta_rad, sigma_l, self.payne_correction)
            f_tp_last = st.f_tp
            HL_last   = st.holdup
            # ΔP_f = f_tp · rho_ns · v² · L / (2·D)  →  v = √(2·D·ΔP / (f_tp·rho_ns·L))
            denom = max(st.f_tp * st.rho_ns * L / D, 1e-30)
            v_m_new = math.sqrt(2.0 * net_driving / denom)
            if abs(v_m_new - v_m) < 1e-5 * max(v_m_new, _MIN_VELOCITY):
                v_m = v_m_new
                break
            v_m = 0.5 * v_m + 0.5 * v_m_new

        return v_m, f_tp_last, HL_last

    def calculate_velocity(
        self,
        pipe_state,
        delta_p: float,
        density: float,
        viscosity: float,
        tolerance: float | None = None,
        colebrook_friction_strategy: str = "transformed",
        friction_factor_method: str = "newton",
        friction_factor_max_iterations: int = 50,
        velocity_loop_method: str = "fixed_point",
        velocity_loop_max_iterations: int = 50,
        velocity_loop_tolerance: float | None = None,
        laminar_turbulent_transition_re: float = 2300.0,
    ) -> float:
        """Return mixture velocity for the given nodal pressure drop.

        The Colebrook-specific keyword arguments (friction_factor_method, etc.)
        are accepted for interface compatibility but are ignored — Beggs-Brill
        uses its own fixed-point iteration.
        """
        D         = pipe_state.component.diameter_m
        L         = pipe_state.component.length_m
        roughness = pipe_state.component.absolute_roughness_m
        h         = pipe_state.component.height_change_m
        theta_rad = math.asin(max(-1.0, min(1.0, h / max(L, 1e-6))))

        C_L, rho_l, rho_g, mu_l, mu_g, sigma_l = self._phase_properties(pipe_state)

        v0 = max(abs(pipe_state.velocity_m_per_s), _MIN_VELOCITY)

        def _net_driving(direction: float):
            """Friction-available pressure drop for flow in ``direction``
            (+1: start → end).  Gravity uses the holdup of that flow direction."""
            st = beggs_brill_gradient(v0, C_L, rho_l, rho_g, mu_l, mu_g,
                                      D, roughness, direction * theta_rad,
                                      sigma_l, self.payne_correction)
            return direction * delta_p - st.dP_dz_gravity * L, st

        # Try the direction the pressure difference points to first.  If it
        # cannot overcome gravity (e.g. a riser whose ΔP is below the static
        # head), the fluid flows the other way, driven by gravity.
        sign = 1.0 if delta_p >= 0.0 else -1.0
        net_driving, st0 = _net_driving(sign)
        if net_driving <= 0.0:
            net_rev, st_rev = _net_driving(-sign)
            if net_rev <= 0.0:
                pipe_state.velocity_m_per_s = 0.0
                pipe_state.friction_factor  = st0.f_tp
                return 0.0
            sign, net_driving, st0 = -sign, net_rev, st_rev

        v_m, f_tp, _HL = self._solve_vm(
            net_driving, C_L, rho_l, rho_g, mu_l, mu_g,
            D, L, roughness, sign * theta_rad, v0, sigma_l,
        )

        pipe_state.velocity_m_per_s = sign * v_m
        pipe_state.friction_factor  = f_tp
        return sign * v_m

    def calculate_coupling(
        self,
        pipe_state,
        density: float,
        viscosity: float,
    ) -> float:
        """Linearised dṁ/dΔP for the SIMPLE pressure-correction matrix.

        From  ΔP = f·ρ·v²·L/(2·D)  with ṁ = ρ·A·v, density cancels:
            d(ṁ)/d(ΔP) = -2·A·D / (f·|v|·L)   (same form as ColebrookPipeCorrelation)
        """
        f = max(abs(pipe_state.friction_factor or 0.02), 1e-8)
        v = max(abs(pipe_state.velocity_m_per_s), 1e-12)
        L = pipe_state.component.length_m
        D = pipe_state.component.diameter_m
        return -2.0 * pipe_state.area_m2 * D / (f * v * L)
