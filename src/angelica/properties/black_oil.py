"""
Black-oil three-phase fluid model.

Correlations
------------
Bubble point / Rs / Bo  : Standing (1947)
Live-oil viscosity      : Beggs & Robinson (1975)
z-factor                : Hall-Yarborough (1974) with Sutton (1985) pseudo-crits
Gas viscosity           : Lee, Gonzalez & Eakin (1966)
Water FVF               : McCain (1990), simplified
Dead-oil properties     : see properties.dead_oil

Sign convention
---------------
All pressures are in Pa, temperatures in °C, densities in kg/m³,
viscosities in Pa·s throughout the Python API.  Internal conversions to
field units (psia, °F, scf/STB) are confined to the correlation functions.
"""
from __future__ import annotations

import functools
import math
from dataclasses import dataclass

from .base import FluidModel

# ── Constants ─────────────────────────────────────────────────────────────────
_P_SC_PA = 101_325.0       # Pa  — standard pressure
_T_SC_K  = 288.706         # K   — standard temperature (60 °F)
_R       = 8.314           # J/(mol·K)
_AIR_DENSITY_SC_KG_M3 = 1.2250  # kg/m³ at (_T_SC_K, _P_SC_PA)

# Conversion factors
_M3M3_TO_SCFSTB = 5.6146   # 1 m³_gas/m³_oil = 5.6146 scf/STB
_PSIA_TO_PA     = 6_894.757
_PA_TO_PSIA     = 1.0 / _PSIA_TO_PA


# ── Unit helpers ──────────────────────────────────────────────────────────────

def _c_to_f(t_c: float) -> float:
    return t_c * 9.0 / 5.0 + 32.0


def _c_to_r(t_c: float) -> float:
    """Celsius to Rankine."""
    return (t_c + 273.15) * 1.8


# ── Pseudo-critical properties (Sutton 1985) ──────────────────────────────────

def _pseudo_critical(gas_gravity: float) -> tuple[float, float]:
    """Return (T_pc [°R], P_pc [psia]) from gas gravity (air = 1)."""
    T_pc = 169.2 + 349.5 * gas_gravity - 74.0 * gas_gravity ** 2
    P_pc = 756.8 - 131.0 * gas_gravity - 3.6  * gas_gravity ** 2
    return T_pc, P_pc


# ── z-factor — Hall-Yarborough (1974) ─────────────────────────────────────────

def z_factor_hall_yarborough(
    pressure_pa: float,
    temperature_c: float,
    gas_gravity: float,
) -> float:
    """Gas compressibility factor via Hall & Yarborough (1974).

    Solves the implicit Hall-Yarborough equation with Newton's method.
    Typical error < 1 % for dry natural gas, 1.05 ≤ T_pr ≤ 3.0.

    Args:
        pressure_pa: Absolute pressure (Pa).
        temperature_c: Temperature (°C).
        gas_gravity: Gas specific gravity (air = 1).
    """
    T_pc, P_pc = _pseudo_critical(gas_gravity)
    T_pr = _c_to_r(temperature_c) / T_pc
    P_pr = pressure_pa * _PA_TO_PSIA / P_pc

    t  = 1.0 / T_pr
    t2 = t * t
    t3 = t * t2
    A  = 0.06125 * t * math.exp(-1.2 * (1.0 - t) ** 2)

    # y is the reduced packing fraction — must stay in (0, 1)
    y = 0.001

    for _ in range(60):
        y  = max(min(y, 0.999), 1e-8)
        y2 = y * y
        y3 = y * y2
        y4 = y * y3
        exp_pow  = 2.18 + 2.82 * t
        exp_term = (90.7*t - 242.2*t2 + 42.4*t3) * y ** exp_pow

        fy = (
            -A * P_pr
            + (y + y2 + y3 - y4) / (1.0 - y) ** 3
            - (14.76*t - 9.76*t2 + 4.58*t3) * y2
            + exp_term
        )
        dfy = (
            (1.0 + 4.0*y + 4.0*y2 - 4.0*y3 + y4) / (1.0 - y) ** 4
            - (29.52*t - 19.52*t2 + 9.16*t3) * y
            + exp_pow * exp_term / y
        )
        if abs(dfy) < 1e-30:
            break
        dy = -fy / dfy
        y  = max(min(y + dy, 0.999), 1e-8)
        if abs(dy) < 1e-10:
            break

    return A * P_pr / y


# ── Bubble point — Standing (1947) ────────────────────────────────────────────

def bubble_point_pa(
    gor_sc_m3_per_m3: float,
    gas_gravity: float,
    api_gravity: float,
    temperature_c: float,
) -> float:
    """Bubble point pressure (Pa) via Standing (1947).

    Args:
        gor_sc_m3_per_m3: Producing GOR at standard conditions (m³/m³).
        gas_gravity: Gas specific gravity (air = 1).
        api_gravity: Stock-tank oil API gravity (°API).
        temperature_c: Temperature (°C).
    """
    if gor_sc_m3_per_m3 <= 0.0:
        return 0.0
    T_F          = _c_to_f(temperature_c)
    Rs_scf_STB   = gor_sc_m3_per_m3 * _M3M3_TO_SCFSTB
    Pb_psia      = 18.2 * (
        (Rs_scf_STB / gas_gravity) ** 0.83
        * 10.0 ** (0.00091 * T_F - 0.0125 * api_gravity)
        - 1.4
    )
    return max(Pb_psia, 0.0) * _PSIA_TO_PA


# ── Solution GOR — Standing (1947) ────────────────────────────────────────────

def solution_gor_m3_per_m3(
    pressure_pa: float,
    temperature_c: float,
    gas_gravity: float,
    api_gravity: float,
    gor_sc_m3_per_m3: float,
) -> float:
    """Solution GOR (m³_gas_sc / m³_oil_sc) via Standing (1947).

    Returns GOR dissolved in oil at (P, T).  Clamped to [0, gor_sc]:
    above the bubble point Rs = gor_sc (all gas dissolved).
    """
    if gor_sc_m3_per_m3 <= 0.0:
        return 0.0
    T_F        = _c_to_f(temperature_c)
    P_psia     = pressure_pa * _PA_TO_PSIA
    Rs_scf_STB = (
        gas_gravity
        * (P_psia / 18.2 + 1.4) ** 1.205
        * 10.0 ** (0.0125 * api_gravity - 0.00091 * T_F)
    )
    Rs_m3_per_m3 = Rs_scf_STB / _M3M3_TO_SCFSTB
    return min(max(Rs_m3_per_m3, 0.0), gor_sc_m3_per_m3)


# ── Oil FVF — Standing (1947) ─────────────────────────────────────────────────

def oil_fvf(
    rs_m3_per_m3: float,
    temperature_c: float,
    gas_gravity: float,
    api_gravity: float,
) -> float:
    """Oil formation volume factor Bo (m³_res / m³_sc) via Standing (1947).

    Args:
        rs_m3_per_m3: Solution GOR at current conditions (m³/m³).
        temperature_c: Temperature (°C).
        gas_gravity: Gas specific gravity (air = 1).
        api_gravity: Stock-tank oil API gravity (°API).
    """
    T_F         = _c_to_f(temperature_c)
    Rs_scf_STB  = rs_m3_per_m3 * _M3M3_TO_SCFSTB
    oil_gravity = 141.5 / (api_gravity + 131.5)
    F = Rs_scf_STB * (gas_gravity / oil_gravity) ** 0.5 + 1.25 * T_F
    return 0.972 + 0.000147 * F ** 1.175


# ── Gas FVF ───────────────────────────────────────────────────────────────────

def gas_fvf(pressure_pa: float, temperature_c: float, z: float) -> float:
    """Gas formation volume factor Bg (m³_res / m³_sc).

    Bg = z * P_sc * T / (T_sc * P)

    Args:
        pressure_pa: Pressure (Pa).
        temperature_c: Temperature (°C).
        z: Gas compressibility factor (–).
    """
    T_K = temperature_c + 273.15
    return z * _P_SC_PA * T_K / (_T_SC_K * pressure_pa)


# ── Water FVF — McCain (1990) ─────────────────────────────────────────────────

def water_fvf(pressure_pa: float, temperature_c: float) -> float:
    """Water formation volume factor Bw (m³_res / m³_sc).

    Simplified McCain (1990) correlation.  Valid for fresh water and
    typical oilfield brine at moderate temperatures and pressures.
    """
    T_F   = _c_to_f(temperature_c)
    P_psia = pressure_pa * _PA_TO_PSIA
    dT    = T_F - 60.0
    Bw    = 1.0 + 1.21e-4 * dT + 1.0e-6 * dT ** 2 - 3.33e-6 * P_psia
    return max(Bw, 0.9)  # physical lower bound


# ── Live-oil viscosity — Beggs & Robinson (1975) ──────────────────────────────

def live_oil_viscosity_pa_s(dead_oil_visc_pa_s: float, rs_m3_per_m3: float) -> float:
    """Saturated (live) oil viscosity via Beggs & Robinson (1975).

    Args:
        dead_oil_visc_pa_s: Dead-oil viscosity at the same temperature (Pa·s).
        rs_m3_per_m3: Solution GOR at current conditions (m³/m³).
    """
    mu_od_cp   = dead_oil_visc_pa_s * 1000.0
    Rs_scf_STB = rs_m3_per_m3 * _M3M3_TO_SCFSTB
    A = 10.715 * (Rs_scf_STB + 100.0) ** (-0.515)
    B = 5.44   * (Rs_scf_STB + 150.0) ** (-0.338)
    mu_o_cp = A * max(mu_od_cp, 1e-6) ** B
    return max(mu_o_cp, 1e-6) * 1e-3  # cP → Pa·s


# ── Gas viscosity — Lee, Gonzalez & Eakin (1966) ─────────────────────────────

def gas_viscosity_pa_s(
    pressure_pa: float,
    temperature_c: float,
    gas_gravity: float,
) -> float:
    """Gas dynamic viscosity via Lee, Gonzalez & Eakin (1966).

    Computes gas density from Hall-Yarborough z-factor internally.

    Args:
        pressure_pa: Pressure (Pa).
        temperature_c: Temperature (°C).
        gas_gravity: Gas specific gravity (air = 1).
    """
    M_g   = gas_gravity * 28.97          # g/mol
    T_R   = _c_to_r(temperature_c)       # Rankine
    T_K   = temperature_c + 273.15
    z     = z_factor_hall_yarborough(pressure_pa, temperature_c, gas_gravity)
    # Gas density in g/cm³
    rho_gcc = pressure_pa * M_g * 1e-3 / (z * _R * T_K) / 1000.0
    K = (9.4 + 0.02 * M_g) * T_R ** 1.5 / (209.0 + 19.0 * M_g + T_R)
    X = 3.5 + 986.0 / T_R + 0.01 * M_g
    Y = 2.4 - 0.2 * X
    mu_cP = K * math.exp(X * max(rho_gcc, 0.0) ** Y) * 1e-4
    return max(mu_cP, 1e-6) * 1e-3  # cP → Pa·s


# ── Water viscosity ───────────────────────────────────────────────────────────

def water_viscosity_pa_s(temperature_c: float) -> float:
    """Water dynamic viscosity (Pa·s) — polynomial fit, valid 0–100 °C.

    Fitted to IAPWS tabulated values; error < 2 % in 0–100 °C range.
    Returns ≈ 1.002e-3 Pa·s at 20 °C (water at standard conditions).
    """
    return math.exp(0.52 - 2.84e-2 * temperature_c + 1.16e-4 * temperature_c ** 2) * 1e-3


# ── Black-oil composition ─────────────────────────────────────────────────────

@dataclass(frozen=True)
class BlackOilComposition:
    """Four-parameter black-oil fluid composition (reservoir characterisation).

    Separates the compositional description of a fluid from the PVT machinery,
    allowing each inlet node in a network to carry its own composition.

    Args:
        api_gravity: Stock-tank oil API gravity (°API).
        gas_gravity: Gas specific gravity relative to air (–).
        gor_sc_m3_per_m3: Producing gas-oil ratio at standard conditions (m³/m³).
        wor_sc_m3_per_m3: Water-oil ratio at standard conditions (m³/m³).
    """
    api_gravity: float
    gas_gravity: float
    gor_sc_m3_per_m3: float
    wor_sc_m3_per_m3: float

    def mix(self, other: "BlackOilComposition", weight_self: float, weight_other: float) -> "BlackOilComposition":
        """Return a mass-weighted mixture of two compositions."""
        w_tot = weight_self + weight_other
        if w_tot <= 0.0:
            return self
        f_s = weight_self  / w_tot
        f_o = weight_other / w_tot
        return BlackOilComposition(
            api_gravity       = f_s * self.api_gravity        + f_o * other.api_gravity,
            gas_gravity       = f_s * self.gas_gravity        + f_o * other.gas_gravity,
            gor_sc_m3_per_m3  = f_s * self.gor_sc_m3_per_m3  + f_o * other.gor_sc_m3_per_m3,
            wor_sc_m3_per_m3  = f_s * self.wor_sc_m3_per_m3  + f_o * other.wor_sc_m3_per_m3,
        )


# ── Standalone PVT function ───────────────────────────────────────────────────

@functools.lru_cache(maxsize=512)
def compute_pvt(
    pressure_pa: float,
    temperature_c: float,
    api_gravity: float,
    gas_gravity: float,
    gor_sc_m3_per_m3: float,
    wor_sc_m3_per_m3: float,
) -> "BlackOilPVTState":
    """Compute the full black-oil PVT state for an explicit fluid composition.

    This is the functional core used by :class:`BlackOilFluid` and by the
    solver's per-pipe composition propagation.  All parameters are explicit so
    the function can be called with any composition, not just the one stored
    in a ``BlackOilFluid`` instance.
    """
    from .dead_oil import (
        dead_oil_density_kg_per_m3,
        dead_oil_specific_heat_j_per_kg_k,
        dead_oil_thermal_conductivity_w_per_m_k,
        dead_oil_viscosity_pa_s,
    )

    P = max(pressure_pa, 1.0)
    T = temperature_c

    Pb = bubble_point_pa(gor_sc_m3_per_m3, gas_gravity, api_gravity, T)
    Rs = solution_gor_m3_per_m3(P, T, gas_gravity, api_gravity, gor_sc_m3_per_m3)
    Bo = oil_fvf(Rs, T, gas_gravity, api_gravity)

    has_gas = gor_sc_m3_per_m3 > 0.0
    if has_gas:
        z  = z_factor_hall_yarborough(P, T, gas_gravity)
        Bg = gas_fvf(P, T, z)
    else:
        z  = 1.0
        Bg = 1.0

    Bw = water_fvf(P, T)

    rho_oil_sc   = dead_oil_density_kg_per_m3(api_gravity)
    rho_gas_sc   = gas_gravity * _AIR_DENSITY_SC_KG_M3
    rho_water_sc = 1_025.0

    rho_oil = (rho_oil_sc + Rs * rho_gas_sc) / Bo
    rho_gas = rho_gas_sc / Bg if has_gas else 0.0
    rho_wtr = rho_water_sc / Bw

    v_oil = Bo
    v_gas = max(0.0, gor_sc_m3_per_m3 - Rs) * Bg if has_gas else 0.0
    v_wtr = wor_sc_m3_per_m3 * Bw
    v_tot = v_oil + v_gas + v_wtr

    alpha_oil = v_oil / v_tot
    alpha_gas = v_gas / v_tot
    alpha_wtr = v_wtr / v_tot

    rho_m = alpha_oil * rho_oil + alpha_gas * rho_gas + alpha_wtr * rho_wtr

    mu_od  = dead_oil_viscosity_pa_s(api_gravity, T)
    mu_oil = live_oil_viscosity_pa_s(mu_od, Rs)
    mu_gas = gas_viscosity_pa_s(P, T, gas_gravity) if has_gas else 0.0
    mu_wtr = water_viscosity_pa_s(T)
    mu_m   = alpha_oil * mu_oil + alpha_gas * mu_gas + alpha_wtr * mu_wtr

    cp_oil = dead_oil_specific_heat_j_per_kg_k(api_gravity, T)
    cp_gas = 2_200.0
    cp_wtr = 4_182.0

    k_oil = dead_oil_thermal_conductivity_w_per_m_k(api_gravity, T)
    k_gas = 0.035
    k_wtr = 0.62

    w_oil = alpha_oil * rho_oil / rho_m
    w_gas = alpha_gas * rho_gas / rho_m if has_gas else 0.0
    w_wtr = alpha_wtr * rho_wtr / rho_m
    cp_m  = w_oil * cp_oil + w_gas * cp_gas + w_wtr * cp_wtr
    k_m   = alpha_oil * k_oil + alpha_gas * k_gas + alpha_wtr * k_wtr

    return BlackOilPVTState(
        pressure_pa=P,
        temperature_c=T,
        bubble_point_pa=Pb,
        rs_m3_per_m3=Rs,
        bo=Bo,
        bg=Bg,
        bw=Bw,
        z=z,
        holdup_oil=alpha_oil,
        holdup_gas=alpha_gas,
        holdup_water=alpha_wtr,
        density_oil_kg_per_m3=rho_oil,
        density_gas_kg_per_m3=rho_gas,
        density_water_kg_per_m3=rho_wtr,
        viscosity_oil_pa_s=mu_oil,
        viscosity_gas_pa_s=mu_gas,
        viscosity_water_pa_s=mu_wtr,
        mixture_density_kg_per_m3=rho_m,
        mixture_viscosity_pa_s=mu_m,
        mixture_specific_heat_j_per_kg_k=cp_m,
        mixture_thermal_conductivity_w_per_m_k=k_m,
    )


# ── PVT state dataclass ───────────────────────────────────────────────────────

@dataclass
class BlackOilPVTState:
    """Complete PVT state for a black-oil mixture at a given (P, T).

    All densities in kg/m³, viscosities in Pa·s, FVFs in m³_res/m³_sc,
    GOR in m³_gas_sc/m³_oil_sc, pressures in Pa.
    """
    pressure_pa: float
    temperature_c: float
    bubble_point_pa: float
    rs_m3_per_m3: float
    bo: float
    bg: float
    bw: float
    z: float
    holdup_oil: float
    holdup_gas: float
    holdup_water: float
    density_oil_kg_per_m3: float
    density_gas_kg_per_m3: float
    density_water_kg_per_m3: float
    viscosity_oil_pa_s: float
    viscosity_gas_pa_s: float
    viscosity_water_pa_s: float
    mixture_density_kg_per_m3: float
    mixture_viscosity_pa_s: float
    mixture_specific_heat_j_per_kg_k: float
    mixture_thermal_conductivity_w_per_m_k: float

    @property
    def undersaturated(self) -> bool:
        """True when P ≥ Pb — no free gas phase."""
        return self.pressure_pa >= self.bubble_point_pa


# ── BlackOilFluid ─────────────────────────────────────────────────────────────

@dataclass(frozen=True)
class BlackOilFluid(FluidModel):
    """Black-oil three-phase fluid model (gas + oil + water).

    Assumes no-slip (homogeneous) flow — all phases travel at the same
    velocity.  Phase fractions are computed from the PVT FVFs and the
    user-specified GOR and WOR at standard conditions.

    Correlations used
    -----------------
    - Bubble point, Rs, Bo : Standing (1947)
    - z-factor             : Hall-Yarborough (1974) + Sutton (1985) pseudo-crits
    - Gas viscosity        : Lee, Gonzalez & Eakin (1966)
    - Live-oil viscosity   : Beggs & Robinson (1975)
    - Dead-oil properties  : see properties.dead_oil
    - Water FVF            : McCain (1990), simplified
    - Water viscosity      : simple exponential fit

    Args:
        api_gravity: Stock-tank oil API gravity (°API).  Typical range 10–58.
        gas_gravity: Gas specific gravity relative to air (–).  Typical 0.55–0.80.
        gor_sc_m3_per_m3: Producing gas-oil ratio at standard conditions
            (m³_gas / m³_oil).  Use 0 for dead-oil (single liquid phase).
        wor_sc_m3_per_m3: Water-oil ratio at standard conditions
            (m³_water / m³_oil).  Use 0 for dry production.
        reference_pressure_pa: Fallback pressure (Pa) when node pressures are
            not yet initialised.  Defaults to 101 325 Pa (1 atm).
        reference_temperature_c: Fallback temperature (°C).  Defaults to 20 °C.
    """

    api_gravity: float
    gas_gravity: float
    gor_sc_m3_per_m3: float
    wor_sc_m3_per_m3: float
    reference_pressure_pa: float = 101_325.0
    reference_temperature_c: float = 20.0

    def __post_init__(self) -> None:
        if self.api_gravity <= -131.5:
            raise ValueError(f"api_gravity must be > -131.5; got {self.api_gravity}")
        if self.gas_gravity <= 0.0 and self.gor_sc_m3_per_m3 > 0.0:
            raise ValueError(f"gas_gravity must be > 0 when GOR > 0; got {self.gas_gravity}")
        if self.gor_sc_m3_per_m3 < 0.0:
            raise ValueError(f"gor_sc_m3_per_m3 must be >= 0; got {self.gor_sc_m3_per_m3}")
        if self.wor_sc_m3_per_m3 < 0.0:
            raise ValueError(f"wor_sc_m3_per_m3 must be >= 0; got {self.wor_sc_m3_per_m3}")

    # ── helpers ───────────────────────────────────────────────────────────────

    def _pressure(self, link_state) -> float:
        start = getattr(link_state, "start_node", None)
        end   = getattr(link_state, "end_node",   None)
        p_s   = getattr(start, "pressure_pa", None) if start is not None else None
        p_e   = getattr(end,   "pressure_pa", None) if end   is not None else None
        if p_s is not None and p_e is not None:
            return 0.5 * (p_s + p_e)
        return p_s if p_s is not None else (p_e if p_e is not None else self.reference_pressure_pa)

    def _temperature(self, link_state) -> float:
        t = getattr(link_state, "temperature_c", None)
        return float(t) if t is not None else self.reference_temperature_c

    # ── PVT ───────────────────────────────────────────────────────────────────

    def pvt(self, pressure_pa: float, temperature_c: float) -> "BlackOilPVTState":
        """Compute the full black-oil PVT state at (pressure_pa, temperature_c)."""
        return compute_pvt(
            pressure_pa,
            temperature_c,
            self.api_gravity,
            self.gas_gravity,
            self.gor_sc_m3_per_m3,
            self.wor_sc_m3_per_m3,
        )

    # ── FluidModel interface ───────────────────────────────────────────────────

    def density_for_link(self, link_state) -> float:
        return self.pvt(
            self._pressure(link_state), self._temperature(link_state)
        ).mixture_density_kg_per_m3

    def viscosity_for_link(self, link_state) -> float:
        return self.pvt(
            self._pressure(link_state), self._temperature(link_state)
        ).mixture_viscosity_pa_s

    def specific_heat_for_link(self, link_state) -> float:
        return self.pvt(
            self._pressure(link_state), self._temperature(link_state)
        ).mixture_specific_heat_j_per_kg_k

    def thermal_conductivity_for_link(self, link_state) -> float:
        return self.pvt(
            self._pressure(link_state), self._temperature(link_state)
        ).mixture_thermal_conductivity_w_per_m_k
