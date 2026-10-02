"""
Beggs-Brill (1973) two-phase pipe-flow pressure-gradient correlation.

Reference
---------
Beggs, H.D., Brill, J.P. (1973). A Study of Two-Phase Flow in Inclined
Pipes. Journal of Petroleum Technology, May, 607-617.

Inclination correction
----------------------
Payne, W.R., Palmer, C.M., Brill, J.P. (1979). Evaluation of Inclined
Pipe, Two-Phase Liquid Holdup and Pressure-Drop Correlations Using
Experimental Data. JPT, Sept., 1198-1208.
(Simplified form: constant correction factors per flow regime, no
surface-tension term required.)

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
from dataclasses import dataclass

from .gravity import GRAVITY_M_PER_S2, elevation_pressure_term
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

def holdup_horizontal(C_L: float, N_Fr: float) -> float:
    """No-slip-corrected liquid holdup H_L(0) for horizontal flow.

    Beggs-Brill (1973) correlations per flow regime.
    H_L(0) is constrained to [C_L, 1] (cannot be less than input fraction).
    """
    if C_L <= 0.0:
        return 0.0
    if C_L >= 1.0:
        return 1.0

    regime, A = flow_regime(C_L, N_Fr)

    if regime == "segregated":
        HL0 = 0.98 * C_L ** 0.4846 / max(N_Fr ** 0.0868, 1e-30)
    elif regime == "intermittent":
        HL0 = 0.845 * C_L ** 0.5351 / max(N_Fr ** 0.0173, 1e-30)
    elif regime == "distributed":
        HL0 = 1.065 * C_L ** 0.5824 / max(N_Fr ** 0.0609, 1e-30)
    else:  # transition
        HL_seg = 0.98  * C_L ** 0.4846 / max(N_Fr ** 0.0868, 1e-30)
        HL_int = 0.845 * C_L ** 0.5351 / max(N_Fr ** 0.0173, 1e-30)
        HL_seg = max(min(HL_seg, 1.0), C_L)
        HL_int = max(min(HL_int, 1.0), C_L)
        return A * HL_seg + (1.0 - A) * HL_int

    return max(min(HL0, 1.0), C_L)


def holdup_inclined(HL0: float, C_L: float, N_Fr: float, theta_rad: float) -> float:
    """Apply Payne et al. (1979) inclination correction to H_L(0).

    Uses simplified constant correction factors (no surface tension required).
    For uphill (θ > 0): holdup is reduced from vertical predictions.
    For downhill (θ < 0): no correction applied (conservative — holdup = H_L(0)).
    """
    if abs(theta_rad) < 1e-6 or C_L <= 0.0 or C_L >= 1.0:
        return HL0

    regime, _ = flow_regime(C_L, N_Fr)

    if theta_rad > 0.0:  # uphill
        if regime == "segregated":
            k = 0.924
        elif regime == "intermittent":
            k = 0.685
        else:  # distributed or transition
            return HL0
        # Sinusoidal angle weighting (Beggs-Brill original form)
        sin_angle = math.sin(1.8 * theta_rad)
        psi = 1.0 + k * (sin_angle - sin_angle ** 3 / 3.0)
        HL = HL0 * psi
    else:  # downhill: conservative — use H_L(0)
        return HL0

    return max(min(HL, 1.0), C_L)


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
    HL0  = holdup_horizontal(C_L, N_Fr)
    HL   = holdup_inclined(HL0, C_L, N_Fr, theta_rad)

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
        regime=regime, holdup=HL, C_L=C_L, N_Fr=N_Fr,
        rho_ns=rho_ns, rho_s=rho_s,
        f_n=f_n, f_tp=f_tp,
        dP_dz_friction=dP_f,
        dP_dz_gravity=dP_g,
        dP_dz_total=dP_total,
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
    - Inclination correction uses Payne et al. simplified factors (no
      surface-tension term); accuracy may be reduced for steep angles.
    - Acceleration pressure drop (small, <1% for most pipelines) is omitted.
    """

    def __init__(self) -> None:
        self._fluid = None

    def set_fluid(self, fluid) -> None:
        """Update the fluid model used to extract phase-fraction data."""
        self._fluid = fluid

    def _phase_properties(self, pipe_state) -> tuple[float, float, float, float, float]:
        """Return (C_L, rho_l, rho_g, mu_l, mu_g) at midpoint conditions."""
        if self._fluid is None or not hasattr(self._fluid, "pvt"):
            return 1.0, 850.0, 1.2, 1e-3, 1.5e-5  # all-liquid fallback

        # Midpoint pressure
        P_s = getattr(getattr(pipe_state, "start_node", None), "pressure_pa", None)
        P_e = getattr(getattr(pipe_state, "end_node",   None), "pressure_pa", None)
        if P_s is not None and P_e is not None:
            P_avg = 0.5 * (float(P_s) + float(P_e))
        elif P_s is not None:
            P_avg = float(P_s)
        elif P_e is not None:
            P_avg = float(P_e)
        else:
            P_avg = getattr(self._fluid, "reference_pressure_pa", 101_325.0)

        T_avg = float(pipe_state.temperature_c) if pipe_state.temperature_c is not None else 20.0

        pvt = self._fluid.pvt(P_avg, T_avg)

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
                                      D, roughness, theta_rad)
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

        C_L, rho_l, rho_g, mu_l, mu_g = self._phase_properties(pipe_state)

        sign = 1.0 if delta_p >= 0.0 else -1.0

        # Gravity contribution at current holdup (iterate once)
        v0 = max(abs(pipe_state.velocity_m_per_s), _MIN_VELOCITY)
        st0 = beggs_brill_gradient(v0, C_L, rho_l, rho_g, mu_l, mu_g,
                                   D, roughness, sign * theta_rad)
        gravity_dp = st0.dP_dz_gravity * L

        net_driving = abs(delta_p) - gravity_dp
        if net_driving <= 0.0:
            pipe_state.velocity_m_per_s = 0.0
            pipe_state.friction_factor  = st0.f_tp
            return 0.0

        v_m, f_tp, _HL = self._solve_vm(
            net_driving, C_L, rho_l, rho_g, mu_l, mu_g,
            D, L, roughness, sign * theta_rad, v0,
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
