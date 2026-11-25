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
from functools import lru_cache

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
    """Bubble point pressure (Pa) via Standing (1947)."""
    if gor_sc_m3_per_m3 <= 0.0:
        return 0.0
    T_F        = _c_to_f(temperature_c)
    Rs_scf_STB = gor_sc_m3_per_m3 * _M3M3_TO_SCFSTB
    Pb_psia    = 18.2 * (
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
    """Solution GOR (m³/m³) via Standing (1947). Clamped to [0, gor_sc]."""
    if gor_sc_m3_per_m3 <= 0.0:
        return 0.0
    T_F        = _c_to_f(temperature_c)
    P_psia     = pressure_pa * _PA_TO_PSIA
    Rs_scf_STB = (
        gas_gravity
        * (P_psia / 18.2 + 1.4) ** 1.205
        * 10.0 ** (0.0125 * api_gravity - 0.00091 * T_F)
    )
    return min(max(Rs_scf_STB / _M3M3_TO_SCFSTB, 0.0), gor_sc_m3_per_m3)


# ── Oil FVF — Standing (1947) ─────────────────────────────────────────────────

def oil_fvf(
    rs_m3_per_m3: float,
    temperature_c: float,
    gas_gravity: float,
    api_gravity: float,
) -> float:
    """Oil FVF Bo (m³_res/m³_sc) via Standing (1947)."""
    T_F         = _c_to_f(temperature_c)
    Rs_scf_STB  = rs_m3_per_m3 * _M3M3_TO_SCFSTB
    oil_gravity = 141.5 / (api_gravity + 131.5)
    F = Rs_scf_STB * (gas_gravity / oil_gravity) ** 0.5 + 1.25 * T_F
    return 0.972 + 0.000147 * F ** 1.175


# ── Bubble point — Vasquez & Beggs (1980) ────────────────────────────────────

def bubble_point_pa_vasquez_beggs(
    gor_sc_m3_per_m3: float,
    gas_gravity: float,
    api_gravity: float,
    temperature_c: float,
) -> float:
    """Bubble point pressure (Pa) via Vasquez & Beggs (1980)."""
    if gor_sc_m3_per_m3 <= 0.0:
        return 0.0
    T_F        = _c_to_f(temperature_c)
    Rs_scf_STB = gor_sc_m3_per_m3 * _M3M3_TO_SCFSTB
    # Gravity correction to separator conditions (100 psia separator)
    gg_norm = gas_gravity * (1.0 + 5.912e-5 * api_gravity * T_F * math.log10(100.0 / 114.7))
    if api_gravity <= 30.0:
        C1, C2, C3 = 0.0362, 1.0937, 25.7240
    else:
        C1, C2, C3 = 0.0178, 1.1870, 23.9310
    Pb_psia = (Rs_scf_STB / (C1 * gg_norm * math.exp(C3 * api_gravity / (T_F + 459.67)))) ** (1.0 / C2)
    return max(Pb_psia, 0.0) * _PSIA_TO_PA


def solution_gor_m3_per_m3_vasquez_beggs(
    pressure_pa: float,
    temperature_c: float,
    gas_gravity: float,
    api_gravity: float,
    gor_sc_m3_per_m3: float,
) -> float:
    """Solution GOR (m³/m³) via Vasquez & Beggs (1980). Clamped to [0, gor_sc]."""
    if gor_sc_m3_per_m3 <= 0.0:
        return 0.0
    T_F    = _c_to_f(temperature_c)
    P_psia = pressure_pa * _PA_TO_PSIA
    T_R    = T_F + 459.67
    gg_norm = gas_gravity * (1.0 + 5.912e-5 * api_gravity * T_F * math.log10(100.0 / 114.7))
    if api_gravity <= 30.0:
        C1, C2, C3 = 0.0362, 1.0937, 25.7240
    else:
        C1, C2, C3 = 0.0178, 1.1870, 23.9310
    Rs_scf_STB = C1 * gg_norm * P_psia ** C2 * math.exp(C3 * api_gravity / T_R)
    return min(max(Rs_scf_STB / _M3M3_TO_SCFSTB, 0.0), gor_sc_m3_per_m3)


def oil_fvf_vasquez_beggs(
    rs_m3_per_m3: float,
    temperature_c: float,
    gas_gravity: float,
    api_gravity: float,
) -> float:
    """Oil FVF Bo (m³_res/m³_sc) via Vasquez & Beggs (1980)."""
    T_F        = _c_to_f(temperature_c)
    Rs_scf_STB = rs_m3_per_m3 * _M3M3_TO_SCFSTB
    if api_gravity <= 30.0:
        C1, C2, C3 = 4.677e-4, 1.751e-5, -1.811e-8
    else:
        C1, C2, C3 = 4.670e-4, 1.100e-5,  1.337e-9
    return 1.0 + C1 * Rs_scf_STB + (T_F - 60.0) * (api_gravity / gas_gravity) * (C2 + C3 * Rs_scf_STB)


# ── Bubble point — Al-Marhoun (1988) ─────────────────────────────────────────

def bubble_point_pa_al_marhoun(
    gor_sc_m3_per_m3: float,
    gas_gravity: float,
    api_gravity: float,
    temperature_c: float,
) -> float:
    """Bubble point pressure (Pa) via Al-Marhoun (1988)."""
    if gor_sc_m3_per_m3 <= 0.0:
        return 0.0
    T_R        = _c_to_r(temperature_c)
    Rs_scf_STB = gor_sc_m3_per_m3 * _M3M3_TO_SCFSTB
    oil_gravity = 141.5 / (api_gravity + 131.5)
    Pb_psia = 5.38088e-3 * Rs_scf_STB ** 0.715082 * gas_gravity ** (-1.87784) * oil_gravity ** 3.1437 * T_R ** 1.32657
    return max(Pb_psia, 0.0) * _PSIA_TO_PA


def solution_gor_m3_per_m3_al_marhoun(
    pressure_pa: float,
    temperature_c: float,
    gas_gravity: float,
    api_gravity: float,
    gor_sc_m3_per_m3: float,
) -> float:
    """Solution GOR (m³/m³) via Al-Marhoun (1988). Clamped to [0, gor_sc]."""
    if gor_sc_m3_per_m3 <= 0.0:
        return 0.0
    T_R         = _c_to_r(temperature_c)
    P_psia      = pressure_pa * _PA_TO_PSIA
    oil_gravity = 141.5 / (api_gravity + 131.5)
    # Invert the Al-Marhoun Pb equation for Rs
    Rs_scf_STB = (
        P_psia / (5.38088e-3 * gas_gravity ** (-1.87784) * oil_gravity ** 3.1437 * T_R ** 1.32657)
    ) ** (1.0 / 0.715082)
    return min(max(Rs_scf_STB / _M3M3_TO_SCFSTB, 0.0), gor_sc_m3_per_m3)


def oil_fvf_al_marhoun(
    rs_m3_per_m3: float,
    temperature_c: float,
    gas_gravity: float,
    api_gravity: float,
) -> float:
    """Oil FVF Bo (m³_res/m³_sc) via Al-Marhoun (1988)."""
    T_R        = _c_to_r(temperature_c)
    Rs_scf_STB = rs_m3_per_m3 * _M3M3_TO_SCFSTB
    oil_gravity = 141.5 / (api_gravity + 131.5)
    F = Rs_scf_STB ** 0.74239 * gas_gravity ** 0.323294 * oil_gravity ** (-1.20204)
    return 0.497069 + 8.62963e-4 * T_R + 1.82594e-3 * F + 3.18099e-6 * F ** 2


# ── Bubble point — Glaso (1980) ───────────────────────────────────────────────

def bubble_point_pa_glaso(
    gor_sc_m3_per_m3: float,
    gas_gravity: float,
    api_gravity: float,
    temperature_c: float,
) -> float:
    """Bubble point pressure (Pa) via Glaso (1980)."""
    if gor_sc_m3_per_m3 <= 0.0:
        return 0.0
    T_F        = _c_to_f(temperature_c)
    Rs_scf_STB = gor_sc_m3_per_m3 * _M3M3_TO_SCFSTB
    X = (Rs_scf_STB / gas_gravity) ** 0.816 * T_F ** 0.172 / api_gravity ** 0.989
    log_X = math.log10(max(X, 1e-30))
    Pb_psia = 10.0 ** (1.7669 + 1.7447 * log_X - 0.30218 * log_X ** 2)
    return max(Pb_psia, 0.0) * _PSIA_TO_PA


def solution_gor_m3_per_m3_glaso(
    pressure_pa: float,
    temperature_c: float,
    gas_gravity: float,
    api_gravity: float,
    gor_sc_m3_per_m3: float,
) -> float:
    """Solution GOR (m³/m³) via Glaso (1980). Clamped to [0, gor_sc]."""
    if gor_sc_m3_per_m3 <= 0.0:
        return 0.0
    T_F    = _c_to_f(temperature_c)
    P_psia = pressure_pa * _PA_TO_PSIA
    # Solve Glaso Pb equation for Rs: Pb = 10^(A + B*logX - C*logX²), X = (Rs/gg)^0.816 * T_F^0.172 / API^0.989
    # Given P_psia, invert for X then for Rs.
    log_P = math.log10(max(P_psia, 1e-6))
    # Quadratic in log_X: 0.30218 logX² - 1.7447 logX + (1.7669 - log_P) = 0
    a_q, b_q, c_q = 0.30218, -1.7447, 1.7669 - log_P
    discriminant = b_q ** 2 - 4.0 * a_q * c_q
    if discriminant < 0.0:
        log_X = -b_q / (2.0 * a_q)
    else:
        log_X = (-b_q - math.sqrt(discriminant)) / (2.0 * a_q)
    X = 10.0 ** log_X
    Rs_scf_STB = gas_gravity * (X * api_gravity ** 0.989 / T_F ** 0.172) ** (1.0 / 0.816)
    return min(max(Rs_scf_STB / _M3M3_TO_SCFSTB, 0.0), gor_sc_m3_per_m3)


def oil_fvf_glaso(
    rs_m3_per_m3: float,
    temperature_c: float,
    gas_gravity: float,
    api_gravity: float,
) -> float:
    """Oil FVF Bo (m³_res/m³_sc) via Glaso (1980)."""
    T_F         = _c_to_f(temperature_c)
    Rs_scf_STB  = rs_m3_per_m3 * _M3M3_TO_SCFSTB
    oil_gravity = 141.5 / (api_gravity + 131.5)
    F = Rs_scf_STB * (gas_gravity / oil_gravity) ** 0.526 + 0.968 * T_F
    log_F = math.log10(max(F, 1e-30))
    Bo_star = 10.0 ** (-6.58511 + 2.91329 * log_F - 0.27683 * log_F ** 2)
    return 1.0 + Bo_star


# ── PVT correlation dispatcher ────────────────────────────────────────────────

_CORRELATIONS: dict[str, dict] = {
    "STANDING": {
        "pb":  bubble_point_pa,
        "rs":  solution_gor_m3_per_m3,
        "bo":  oil_fvf,
    },
    "VASQUEZ_BEGGS": {
        "pb":  bubble_point_pa_vasquez_beggs,
        "rs":  solution_gor_m3_per_m3_vasquez_beggs,
        "bo":  oil_fvf_vasquez_beggs,
    },
    "AL_MARHOUN": {
        "pb":  bubble_point_pa_al_marhoun,
        "rs":  solution_gor_m3_per_m3_al_marhoun,
        "bo":  oil_fvf_al_marhoun,
    },
    "GLASO": {
        "pb":  bubble_point_pa_glaso,
        "rs":  solution_gor_m3_per_m3_glaso,
        "bo":  oil_fvf_glaso,
    },
}

SUPPORTED_BLACK_OIL_CORRELATIONS: dict[str, str] = {
    "STANDING":      "Standing (1947) — default, valid worldwide",
    "VASQUEZ_BEGGS": "Vasquez & Beggs (1980) — good for API 15–60",
    "AL_MARHOUN":    "Al-Marhoun (1988) — Middle-East oils",
    "GLASO":         "Glaso (1980) — North Sea oils",
}


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


@dataclass(frozen=True)
class GasStream:
    """A stream of gas only (no stock-tank oil or water), e.g. lift gas, for
    :func:`mix_black_oil_streams`."""
    gas_gravity: float


_GAS_ONLY_GOR = 1.0e6   # GOR given to a node that receives gas but no oil


def mix_black_oil_streams(contributions, default: BlackOilComposition) -> BlackOilComposition:
    """Black-oil composition of the mixture of ``(mass_flow, composition)``
    streams, conserving standard volumes: the oil, gas and water standard
    volumes add; GOR and WOR are per mixed oil volume, the gas gravity is
    gas-volume weighted and the oil specific gravity oil-volume weighted.
    ``composition`` is a :class:`BlackOilComposition` or a :class:`GasStream`.

    (Before v1.25.0 the four parameters were averaged by mass, which does not
    conserve the standard volumes and cannot add a gas-only stream.)
    A node receiving gas but no oil gets a GOR of 1e6 and ``default``'s API.
    """
    from .dead_oil import dead_oil_density_kg_per_m3

    q_o = q_g = q_w = g_sum = sg_sum = 0.0
    for m, c in contributions:
        if m <= 0.0:
            continue
        if isinstance(c, GasStream):
            v_g = m / (c.gas_gravity * _AIR_DENSITY_SC_KG_M3)
            q_g += v_g
            g_sum += v_g * c.gas_gravity
            continue
        rho_o = dead_oil_density_kg_per_m3(c.api_gravity)
        m_per_oil = rho_o + c.gor_sc_m3_per_m3 * c.gas_gravity * _AIR_DENSITY_SC_KG_M3 + c.wor_sc_m3_per_m3 * 1_025.0
        v_o = m / m_per_oil
        q_o += v_o
        q_g += v_o * c.gor_sc_m3_per_m3
        q_w += v_o * c.wor_sc_m3_per_m3
        g_sum += v_o * c.gor_sc_m3_per_m3 * c.gas_gravity
        sg_sum += v_o * rho_o / 999.064          # as in dead_oil_density_kg_per_m3
    if q_o <= 0.0 and q_g <= 0.0:
        return default
    gg = g_sum / q_g if q_g > 0.0 else default.gas_gravity
    if q_o <= 0.0:
        return BlackOilComposition(api_gravity=default.api_gravity, gas_gravity=gg,
                                   gor_sc_m3_per_m3=_GAS_ONLY_GOR, wor_sc_m3_per_m3=0.0)
    sg = sg_sum / q_o
    return BlackOilComposition(api_gravity=141.5 / sg - 131.5, gas_gravity=gg,
                               gor_sc_m3_per_m3=q_g / q_o, wor_sc_m3_per_m3=q_w / q_o)


# ── Standalone PVT function ───────────────────────────────────────────────────

@functools.lru_cache(maxsize=512)
def compute_pvt(
    pressure_pa: float,
    temperature_c: float,
    api_gravity: float,
    gas_gravity: float,
    gor_sc_m3_per_m3: float,
    wor_sc_m3_per_m3: float,
    correlation: str = "STANDING",
    gas_specific_heat_j_per_kg_k: float = 2_200.0,
    gas_thermal_conductivity_w_per_mk: float = 0.035,
) -> "BlackOilPVTState":
    """Compute the full black-oil PVT state for an explicit fluid composition.

    All parameters are explicit so the function can be called with any
    composition, not just the one stored in a ``BlackOilFluid`` instance.
    The ``correlation`` key selects which Pb/Rs/Bo correlation set to use
    (see :data:`SUPPORTED_BLACK_OIL_CORRELATIONS`).
    """
    from .dead_oil import (
        dead_oil_density_kg_per_m3,
        dead_oil_specific_heat_j_per_kg_k,
        dead_oil_thermal_conductivity_w_per_m_k,
        dead_oil_viscosity_pa_s,
    )

    corr = _CORRELATIONS.get(correlation, _CORRELATIONS["STANDING"])

    P = max(pressure_pa, 1.0)
    # Below T_MIN_C the dead-oil viscosity correlation has no value (log of
    # T in °F ≤ 0): properties are taken at T_MIN_C, so enthalpy continues
    # linearly with that cp.  The black-oil solver reports a solution with a
    # temperature below T_MIN_C as not converged (outside the correlations).
    # Above T_MAX_C nothing is physical (it is reached only by the energy
    # solve of an unsettled hydraulic pass); the water viscosity's
    # exp(… + 1.16e-4·T²) would overflow near 2 500 °C.
    T = min(max(temperature_c, T_MIN_C), T_MAX_C)

    Pb = corr["pb"](gor_sc_m3_per_m3, gas_gravity, api_gravity, T)
    Rs = corr["rs"](P, T, gas_gravity, api_gravity, gor_sc_m3_per_m3)
    Bo = corr["bo"](Rs, T, gas_gravity, api_gravity)

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
    cp_gas = gas_specific_heat_j_per_kg_k
    cp_wtr = 4_182.0

    k_oil = dead_oil_thermal_conductivity_w_per_m_k(api_gravity, T)
    k_gas = gas_thermal_conductivity_w_per_mk
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


# ── Enthalpy (energy balance: Joule-Thomson, friction heating) ────────────────
#
# h(P, T) = ∫_0^T cp(P0, T') dT' + ∫_{P0}^{P} (∂h/∂P)_T dP',   P0 = 1 atm,
# (∂h/∂P)_T = Σ_phases w_i · (v_i − T·(∂v_i/∂T)_P)              (T in K)
#
# w_i are the in-situ mass fractions of oil, free gas and water, v_i = 1/ρ_i.
# For the free gas v = Z·R·T/(P·M), so v − T·∂v/∂T = −R·T²/(P·M)·∂Z/∂T: a real
# gas cools on expansion.  For the liquids v dominates: they warm (friction
# heating).  The heat of solution of the dissolved gas is neglected — gas
# released from the oil only changes the phase split.  h is zero at 1 atm and
# 0 °C, which is also the reference of the global energy balance.

_GL4 = ((-0.8611363115940526, 0.3478548451374538), (-0.3399810435848563, 0.6521451548625461),
        (0.3399810435848563, 0.6521451548625461), (0.8611363115940526, 0.3478548451374538))
_GL8 = ((-0.9602898564975363, 0.1012285362903763), (-0.7966664774136267, 0.2223810344533745),
        (-0.5255324099163290, 0.3137066458778873), (-0.1834346424956498, 0.3626837833783620),
        (0.1834346424956498, 0.3626837833783620), (0.5255324099163290, 0.3137066458778873),
        (0.7966664774136267, 0.2223810344533745), (0.9602898564975363, 0.1012285362903763))
_P0_PA = 101_325.0


def _dh_dp(params: tuple, P: float, T_c: float) -> float:
    """(∂h/∂P)_T of the black-oil mixture, m³/kg (= J/(kg·Pa))."""
    dT = 0.5
    s0 = compute_pvt(P, T_c, *params)
    s_lo = compute_pvt(P, T_c - dT, *params)
    s_hi = compute_pvt(P, T_c + dT, *params)
    T_k = T_c + 273.15
    total = 0.0
    for a, rho in (("oil", "density_oil_kg_per_m3"), ("gas", "density_gas_kg_per_m3"),
                   ("water", "density_water_kg_per_m3")):
        alpha = getattr(s0, f"holdup_{a}")
        r0, rl, rh = getattr(s0, rho), getattr(s_lo, rho), getattr(s_hi, rho)
        if alpha <= 0.0 or r0 <= 0.0 or rl <= 0.0 or rh <= 0.0:
            continue
        w = alpha * r0 / s0.mixture_density_kg_per_m3
        dv_dT = (1.0 / rh - 1.0 / rl) / (2.0 * dT)
        total += w * (1.0 / r0 - T_k * dv_dT)
    return total


@lru_cache(maxsize=16384)
def black_oil_enthalpy(params: tuple, pressure_pa: float, temperature_c: float) -> float:
    """Specific enthalpy (J/kg) of a black-oil mixture; ``params`` are the
    trailing arguments of :func:`compute_pvt` (api, gas gravity, GOR, WOR,
    correlation, gas cp, gas k)."""
    T = temperature_c
    h_T = 0.5 * T * sum(w * compute_pvt(_P0_PA, 0.5 * T * (1.0 + x), *params).mixture_specific_heat_j_per_kg_k
                        for x, w in _GL4)
    P = max(pressure_pa, 1.0)
    half = 0.5 * (P - _P0_PA)
    mid = 0.5 * (P + _P0_PA)
    h_P = half * sum(w * _dh_dp(params, mid + half * x, T) for x, w in _GL8)
    return h_T + h_P


# Lowest temperature the black-oil correlations take: the stated lower limit
# of the Beggs-Robinson dead-oil viscosity correlation (−6.7 °C, 20 °F, as
# given in ``angelica.properties.dead_oil``).  Its formula 10^(10^z) − 1
# overflows a few °F above 0 °F for heavy oils, so properties are not
# evaluated below this limit.
T_MIN_C = -6.7
T_MAX_C = 500.0


def temperature_from_enthalpy(h_func, target: float, T_guess: float) -> float:
    """Invert a monotonic h(T) (secant steps safeguarded by bisection)."""
    lo, hi = T_guess - 5.0, T_guess + 5.0
    f_lo, f_hi = h_func(lo) - target, h_func(hi) - target
    for _ in range(40):
        if f_lo * f_hi <= 0.0:
            break
        if abs(f_lo) < abs(f_hi):
            lo -= 2.0 * (hi - lo)
            f_lo = h_func(lo) - target
        else:
            hi += 2.0 * (hi - lo)
            f_hi = h_func(hi) - target
    else:
        return T_guess
    # Regula falsi with the Illinois modification: plain regula falsi keeps
    # one end of the bracket fixed and creeps, and h is evaluated at T
    # rounded to 1e-4 K (steps of ~0.3 J/kg), so |f| < 1e-3 J/kg may never
    # be met.  Until v1.27.0 the loop then returned the bracket's midpoint —
    # up to 2.5 K from the root, alternating between iterations.
    T = T_prev = 0.5 * (lo + hi)
    side = 0
    for _ in range(80):
        T = hi - f_hi * (hi - lo) / (f_hi - f_lo) if f_hi != f_lo else 0.5 * (lo + hi)
        if not (lo < T < hi):
            T = 0.5 * (lo + hi)
        f = h_func(T) - target
        if abs(f) < 1e-3 or hi - lo < 1e-5 or abs(T - T_prev) < 1e-6:
            return T
        T_prev = T
        if f_lo * f <= 0.0:
            hi, f_hi = T, f
            if side == -1:
                f_lo *= 0.5
            side = -1
        else:
            lo, f_lo = T, f
            if side == 1:
                f_hi *= 0.5
            side = 1
    return T


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
    correlation: str = "STANDING"
    gas_specific_heat_j_per_kg_k: float = 2_200.0
    gas_thermal_conductivity_w_per_mk: float = 0.035

    def __post_init__(self) -> None:
        if self.api_gravity <= -131.5:
            raise ValueError(f"api_gravity must be > -131.5; got {self.api_gravity}")
        if self.gas_gravity <= 0.0 and self.gor_sc_m3_per_m3 > 0.0:
            raise ValueError(f"gas_gravity must be > 0 when GOR > 0; got {self.gas_gravity}")
        if self.gor_sc_m3_per_m3 < 0.0:
            raise ValueError(f"gor_sc_m3_per_m3 must be >= 0; got {self.gor_sc_m3_per_m3}")
        if self.wor_sc_m3_per_m3 < 0.0:
            raise ValueError(f"wor_sc_m3_per_m3 must be >= 0; got {self.wor_sc_m3_per_m3}")
        if self.correlation not in SUPPORTED_BLACK_OIL_CORRELATIONS:
            raise ValueError(
                f"Unknown black-oil correlation '{self.correlation}'. "
                f"Valid options: {list(SUPPORTED_BLACK_OIL_CORRELATIONS)}"
            )

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
            self.correlation,
            self.gas_specific_heat_j_per_kg_k,
            self.gas_thermal_conductivity_w_per_mk,
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

    # ── enthalpy (see black_oil_enthalpy) ─────────────────────────────────────

    def _pvt_params(self) -> tuple:
        return (self.api_gravity, self.gas_gravity, self.gor_sc_m3_per_m3, self.wor_sc_m3_per_m3,
                self.correlation, self.gas_specific_heat_j_per_kg_k, self.gas_thermal_conductivity_w_per_mk)

    def enthalpy_for_link(self, link_state, pressure_pa: float, temperature_c: float) -> float:
        """Specific enthalpy (J/kg) at (P, T); zero at 1 atm and 0 °C."""
        return black_oil_enthalpy(self._pvt_params(), float(round(pressure_pa / 10.0) * 10.0),
                                  round(temperature_c, 4))

    def reference_enthalpy_for_link(self, link_state) -> float:
        return 0.0

    def specific_heat_at(self, link_state, pressure_pa: float, temperature_c: float) -> float:
        return self.pvt(pressure_pa, temperature_c).mixture_specific_heat_j_per_kg_k

    def temperature_from_enthalpy_for_link(
        self, link_state, enthalpy_j_per_kg: float, pressure_pa: float, T_guess_c: float = 20.0,
    ) -> float:
        return temperature_from_enthalpy(
            lambda T: self.enthalpy_for_link(link_state, pressure_pa, T), enthalpy_j_per_kg, T_guess_c)
