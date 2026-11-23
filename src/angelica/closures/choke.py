"""
Two-phase choke flow: Sachdeva et al. (1986) critical / subcritical model.

A production choke throttles a gas-liquid mixture through a short restriction.
When the downstream-to-upstream pressure ratio y = P2/P1 falls below a
critical value y_c, the throat velocity reaches the mixture sonic velocity and
the mass flow rate no longer depends on the downstream pressure ("critical"
or "choked" flow).  Above y_c the flow is subcritical and responds to P2.

Model (SI units, g_c = 1)
-------------------------
Upstream state (subscript 1): pressure P1, gas mass fraction x1 (frozen
through the restriction — no interphase mass transfer), gas specific volume
v_G1 = 1/ρ_G1, liquid specific volume v_L = 1/ρ_L (incompressible).  The gas
expands isentropically: v_G2 = v_G1 · y^(-1/k).

Critical pressure ratio y_c solves y = f(y):

    f(y) = [ (k/(k-1) + (1-x1)·v_L·(1-y) / (x1·v_G1))
             / (k/(k-1) + n/2 + n·R2 + (n/2)·R2²) ] ^ (k/(k-1))

    R2 = (1-x1)·v_L / (x1·v_G2),
    n  = 1 + x1·(Cp - Cv) / (x1·Cv + (1-x1)·C_L)    (polytropic exponent)

Mass flux through the throat, evaluated at y = max(y_c, P2/P1):

    G = Cd · ρ_m2 · √( 2·P1·[ (1-x1)·(1-y)·v_L + x1·k/(k-1)·(v_G1 - y·v_G2) ] )

    1/ρ_m2 = x1·v_G2 + (1-x1)·v_L

Limits (verified in tests):
  - x1 → 0 (liquid):  y_c → 0, G → Cd·√(2·ρ_L·ΔP)       (orifice equation)
  - x1 = 1 (gas):     y_c = (2/(k+1))^(k/(k-1))            (isentropic nozzle)

Reference
---------
Sachdeva, R., Schmidt, Z., Brill, J.P. and Blais, R.M. (1986). Two-Phase Flow
Through Chokes. SPE 15657.  Also: Brill, J.P. and Mukherjee, H. (1999).
Multiphase Flow in Wells, SPE Monograph Vol. 17, §5.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from .minor_losses import MinorLossModel

_R_UNIVERSAL = 8_314.462618  # J/(kmol·K)
_M_AIR = 28.9647             # kg/kmol
_X_LIQUID = 1e-9             # gas mass fraction below which flow is treated as liquid


# ── Core correlation ────────────────────────────────────────────────────────────

def _f_ratio(y: float, x1: float, v_g1: float, v_l: float, k: float, n: float) -> float:
    """Right-hand side f(y) of the Sachdeva critical-ratio equation."""
    kk = k / (k - 1.0)
    v_g2 = v_g1 * y ** (-1.0 / k)
    r1 = (1.0 - x1) * v_l / (x1 * v_g1)
    r2 = (1.0 - x1) * v_l / (x1 * v_g2)
    num = kk + r1 * (1.0 - y)
    den = kk + 0.5 * n + n * r2 + 0.5 * n * r2 * r2
    return (num / den) ** kk


def sachdeva_critical_ratio(
    x1: float, rho_g1: float, rho_l: float, k: float, n: float,
    tol: float = 1e-12, max_iter: int = 200,
) -> float:
    """Critical pressure ratio y_c = P2/P1 at which the choke flow becomes sonic.

    Solved by bisection on g(y) = y - f(y), which changes sign on (0, 1):
    g(0⁺) = -f(0⁺) < 0 and g(1) = 1 - f(1) > 0.  Returns 0 for pure liquid.
    """
    if x1 <= _X_LIQUID:
        return 0.0
    v_g1, v_l = 1.0 / rho_g1, 1.0 / rho_l
    lo, hi = 1e-12, 1.0
    for _ in range(max_iter):
        mid = 0.5 * (lo + hi)
        if mid - _f_ratio(mid, x1, v_g1, v_l, k, n) > 0.0:
            hi = mid
        else:
            lo = mid
        if hi - lo < tol:
            break
    return 0.5 * (lo + hi)


def sachdeva_mass_flux(
    p1_pa: float, y: float, x1: float, rho_g1: float, rho_l: float, k: float, cd: float,
) -> float:
    """Throat mass flux G (kg/m²/s) at pressure ratio y (use max(y_c, P2/P1))."""
    y = min(max(y, 0.0), 1.0)
    v_l = 1.0 / rho_l
    if x1 <= _X_LIQUID:
        return cd * math.sqrt(max(2.0 * rho_l * p1_pa * (1.0 - y), 0.0))
    v_g1 = 1.0 / rho_g1
    v_g2 = v_g1 * y ** (-1.0 / k) if y > 0.0 else math.inf
    rho_m2 = 1.0 / (x1 * v_g2 + (1.0 - x1) * v_l)
    work = (1.0 - x1) * (1.0 - y) * v_l + x1 * k / (k - 1.0) * v_g1 * (1.0 - y ** ((k - 1.0) / k))
    return cd * rho_m2 * math.sqrt(max(2.0 * p1_pa * work, 0.0))


@dataclass(frozen=True)
class ChokeOperatingPoint:
    """Result of a choke calculation at given upstream/downstream pressures."""

    mass_flow_kg_per_s: float
    mass_flux_kg_per_m2s: float
    critical_ratio: float        # y_c
    pressure_ratio: float        # actual P2/P1
    critical: bool               # True when P2/P1 <= y_c (flow is choked)
    gas_mass_fraction: float     # x1 at upstream conditions
    heat_capacity_ratio: float   # k = Cp/Cv of the gas


def choke_flow(
    p1_pa: float, p2_pa: float, orifice_area_m2: float, cd: float,
    x1: float, rho_g1: float, rho_l: float, k: float, n: float,
) -> ChokeOperatingPoint:
    """Mass flow through a choke for upstream P1 ≥ downstream P2."""
    y_c = sachdeva_critical_ratio(x1, rho_g1, rho_l, k, n)
    y_act = p2_pa / p1_pa if p1_pa > 0.0 else 1.0
    y_eff = max(y_c, y_act)
    G = sachdeva_mass_flux(p1_pa, y_eff, x1, rho_g1, rho_l, k, cd)
    return ChokeOperatingPoint(
        mass_flow_kg_per_s=G * orifice_area_m2,
        mass_flux_kg_per_m2s=G,
        critical_ratio=y_c,
        pressure_ratio=y_act,
        critical=x1 > _X_LIQUID and y_act <= y_c,
        gas_mass_fraction=x1,
        heat_capacity_ratio=k,
    )


# ── Black-oil property extraction ──────────────────────────────────────────────

def black_oil_choke_inputs(fluid, p1_pa: float, temperature_c: float) -> tuple[float, float, float, float, float]:
    """Return (x1, ρ_G1, ρ_L, k, n) for a BlackOilFluid at upstream conditions.

    k uses ideal-gas Cv = Cp - R/M with M = 28.97·γ_g; the liquid heat
    capacity is recovered from the mass-weighted mixture value.
    """
    pvt = fluid.pvt(p1_pa, temperature_c)
    alpha_l = pvt.holdup_oil + pvt.holdup_water
    rho_l = (
        (pvt.holdup_oil * pvt.density_oil_kg_per_m3 + pvt.holdup_water * pvt.density_water_kg_per_m3)
        / alpha_l if alpha_l > 0.0 else pvt.density_oil_kg_per_m3
    )
    rho_g = pvt.density_gas_kg_per_m3
    x1 = pvt.holdup_gas * rho_g / pvt.mixture_density_kg_per_m3 if pvt.holdup_gas > 0.0 else 0.0

    cp_g = float(getattr(fluid, "gas_specific_heat_j_per_kg_k", 2_200.0))
    gas_gravity = float(getattr(fluid, "gas_gravity", 0.65)) or 0.65
    cv_g = max(cp_g - _R_UNIVERSAL / (_M_AIR * gas_gravity), 0.5 * cp_g)
    k = cp_g / cv_g
    if x1 < 1.0 - 1e-12:
        c_l = max((pvt.mixture_specific_heat_j_per_kg_k - x1 * cp_g) / (1.0 - x1), 1.0)
    else:
        c_l = 2_000.0
    n = 1.0 + x1 * (cp_g - cv_g) / (x1 * cv_g + (1.0 - x1) * c_l)
    return x1, rho_g if rho_g > 0.0 else 1.0, rho_l, k, n


def resolve_link_fluid(fluid, link_state):
    """Return the per-link fluid (handles the solver's per-pipe wrapper)."""
    if fluid is not None and hasattr(fluid, "fluid_for_link"):
        return fluid.fluid_for_link(link_state)
    return fluid


def choke_operating_point(choke_state, fluid) -> ChokeOperatingPoint | None:
    """Evaluate a solved choke (a FittingState wrapping a Choke).

    Returns ``None`` when the fluid does not expose black-oil PVT, in which
    case the solver uses the incompressible orifice equation.
    """
    link_fluid = resolve_link_fluid(fluid, choke_state)
    if link_fluid is None or not hasattr(link_fluid, "pvt"):
        return None
    p_s = float(choke_state.start_node.pressure_pa)
    p_e = float(choke_state.end_node.pressure_pa)
    if p_s >= p_e:
        p1, p2, up_node = p_s, p_e, choke_state.start_node
    else:
        p1, p2, up_node = p_e, p_s, choke_state.end_node
    T = getattr(up_node, "temperature_c", None)
    T = float(T) if T is not None else float(getattr(link_fluid, "reference_temperature_c", 20.0))
    x1, rho_g, rho_l, k, n = black_oil_choke_inputs(link_fluid, p1, T)
    comp = choke_state.component
    area = math.pi * comp.orifice_diameter_m ** 2 / 4.0
    return choke_flow(p1, p2, area, comp.discharge_coefficient, x1, rho_g, rho_l, k, n)


# ── Solver correlation ─────────────────────────────────────────────────────────

class ChokeFlowModel(MinorLossModel):
    """Device model of the black-oil solver.

    - ``Choke`` links: Sachdeva critical / subcritical two-phase flow;
    - ``WellInflow`` links: IPR (``angelica.closures.well_inflow``);
    - ``Well`` links: IPR + production tubing, stable branch of the
      deliverability curve (``angelica.closures.well_model``);
    (wells take a black-oil or a compositional fluid);
    - everything else (ordinary fittings, check valves, orifice-model chokes,
      fluids without black-oil PVT): the K-factor equation of
      ``MinorLossModel``.

    ``SteadyBlackOilSolver`` uses this model by default and calls
    ``set_fluid()`` before every hydraulic solve.
    """

    def __init__(self) -> None:
        self._fluid = None
        # Set by the solver during its turbulent (physical) hydraulic phase:
        # producing wells may then die at the maximum of their curve.
        self.well_latch_active = False

    def set_fluid(self, fluid) -> None:
        self._fluid = fluid

    def _uses_sachdeva(self, link_state) -> bool:
        comp = link_state.component
        if getattr(comp, "flow_model", None) != "sachdeva":
            return False
        link_fluid = resolve_link_fluid(self._fluid, link_state)
        return link_fluid is not None and hasattr(link_fluid, "pvt")

    def _well_fluid(self, link_state, flag: str = "is_well_inflow"):
        if not getattr(link_state.component, flag, False):
            return None
        from angelica.properties.compositional_fluid import CompositionalFluid

        link_fluid = resolve_link_fluid(self._fluid, link_state)
        if link_fluid is None or not (hasattr(link_fluid, "pvt") or isinstance(link_fluid, CompositionalFluid)):
            raise ValueError("WellInflow and Well links need a black-oil or compositional fluid.")
        return link_fluid

    def calculate_velocity(self, link_state, delta_p, density, viscosity, tolerance=None) -> float:
        well_model_fluid = self._well_fluid(link_state, "is_well_model")
        if well_model_fluid is not None:
            from .well_model import well_link_curve
            curve = well_link_curve(link_state, well_model_fluid)
            p_wh = float(link_state.end_node.pressure_pa)
            # A producing well dies when the wellhead reaches the curve's
            # maximum — only in the turbulent (physical) phase; the laminar
            # warm-start pressures are not physical.  Starting again is
            # decided by the solver between outer iterations.
            if (self.well_latch_active and getattr(link_state, "well_latch_enabled", False)
                    and curve.dies_at(p_wh)):
                link_state.well_dead = True
            m_dot, _ = curve.mass_flow(p_wh, getattr(link_state, "well_dead", False))
            return m_dot / (density * link_state.area_m2)
        well_fluid = self._well_fluid(link_state)
        if well_fluid is not None:
            from .well_inflow import well_inflow_point
            m_dot = well_inflow_point(link_state, well_fluid).mass_flow_kg_per_s
            return m_dot / (density * link_state.area_m2)
        if not self._uses_sachdeva(link_state) or delta_p == 0.0:
            return super().calculate_velocity(link_state, delta_p, density, viscosity, tolerance)
        op = choke_operating_point(link_state, self._fluid)
        assert op is not None
        v = op.mass_flow_kg_per_s / (density * link_state.area_m2)
        return v if delta_p > 0.0 else -v

    def calculate_coupling(self, link_state, density, viscosity) -> float:
        well_model_fluid = self._well_fluid(link_state, "is_well_model")
        if well_model_fluid is not None:
            # ΔP = P_res − P_wh with P_res fixed:  dṁ/dΔP = −dṁ/dP_wh  (> 0 on the
            # stable branch).  A dead well gets a near-zero coupling.
            from .well_model import well_link_curve
            curve = well_link_curve(link_state, well_model_fluid)
            _, dm_dpwh = curve.mass_flow(float(link_state.end_node.pressure_pa),
                                         getattr(link_state, "well_dead", False))
            return min(dm_dpwh, -1e-12)
        well_fluid = self._well_fluid(link_state)
        if well_fluid is not None:
            # Tangent dṁ/dΔP of the IPR (exact for the linear part).
            from .well_inflow import well_inflow_coupling
            return well_inflow_coupling(link_state, well_fluid)
        if not self._uses_sachdeva(link_state):
            return super().calculate_coupling(link_state, density, viscosity)
        delta_p = abs(float(link_state.start_node.pressure_pa) - float(link_state.end_node.pressure_pa))
        if delta_p < 1.0:
            return super().calculate_coupling(link_state, density, viscosity)
        # Secant coupling ṁ/ΔP — same convention as the K-factor and pipe models.
        m_dot = density * abs(link_state.velocity_m_per_s) * link_state.area_m2
        return -max(m_dot, 1e-12) / delta_p
