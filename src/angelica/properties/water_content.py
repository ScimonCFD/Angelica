"""
Water content of a gas and its water dew point (roadmap item 2b.1).

The model is Angelica's free-water model: ideal gas and Raoult's law for
the water vapour, y_w·P = P_sat(T) at saturation, with the water saturation
pressure of ``free_water.water_psat_pa`` (Wagner / IAPWS-IF97 form).  It
neglects the non-ideal enhancement of the water content at pressure
(fugacity and Poynting corrections), so it underestimates the water a gas
holds at high pressure — not quantified here (no measured water-content
data were read; see ``tests/test_water_content.py``).

Water content inputs:
  - mole fraction (or ppm by mole);
  - mg per standard m³ of gas (60 °F, 1 atm, ideal gas: 42.21 mol/Sm³);
  - saturation fraction s at (T, P): y_w = s·P_sat(T)/P.

The water dew point at P is the temperature at which the gas is saturated,
P_sat(T_dew) = y_w·P.  Below 0 °C the result is the dew point over
supercooled liquid water (P_sat extrapolated below its 273 K validity);
the frost point (over ice) is not modelled.  Dew points outside
−50 … 373.85 °C raise ValueError instead of being clamped.
"""
from __future__ import annotations

import math

from .free_water import water_psat_pa

M_WATER = 0.018015                     # kg/mol
_P_SC, _T_SC = 101_325.0, 288.706      # 60 °F, 1 atm (Angelica's standard conditions)
_R = 8.314462618


def mol_per_sm3() -> float:
    """Moles of ideal gas per standard m³ (60 °F, 1 atm)."""
    return _P_SC / (_R * _T_SC)


def saturated_water_mole_fraction(temperature_c: float, pressure_pa: float) -> float:
    """Water mole fraction of a gas saturated with water at (T, P)."""
    if not (math.isfinite(pressure_pa) and pressure_pa > 0.0 and math.isfinite(temperature_c)):
        raise ValueError(f"need a finite temperature and a pressure > 0 (got {temperature_c} °C, {pressure_pa} Pa)")
    return min(water_psat_pa(temperature_c + 273.15) / pressure_pa, 1.0)


def water_mole_fraction(*, ppm_mol: float | None = None, mg_per_sm3: float | None = None,
                        saturation: float | None = None, temperature_c: float | None = None,
                        pressure_pa: float | None = None) -> float:
    """Water mole fraction of a gas from one of: ppm by mole, mg per standard
    m³, or a saturation fraction at (``temperature_c``, ``pressure_pa``)."""
    given = [v for v in (ppm_mol, mg_per_sm3, saturation) if v is not None]
    if len(given) != 1:
        raise ValueError("give exactly one of ppm_mol, mg_per_sm3, saturation")
    if not (math.isfinite(given[0]) and given[0] >= 0.0):
        raise ValueError(f"water content must be finite and ≥ 0 (got {given[0]})")
    if ppm_mol is not None:
        return _fraction(ppm_mol * 1e-6)
    if mg_per_sm3 is not None:
        return _fraction(mg_per_sm3 * 1e-6 / M_WATER / mol_per_sm3())
    if saturation is None or temperature_c is None or pressure_pa is None:
        raise ValueError("a saturation fraction needs temperature_c and pressure_pa")
    if saturation > 1.0:
        raise ValueError(f"saturation fraction must be ≤ 1 (got {saturation})")
    return saturation * saturated_water_mole_fraction(temperature_c, pressure_pa)


def _fraction(y: float) -> float:
    if y > 1.0:
        raise ValueError(f"water content gives a mole fraction above 1 ({y:.4g})")
    return y


def water_dew_point_c(water_mole_fraction_gas: float, pressure_pa: float) -> float | None:
    """Water dew point (°C) of a gas at ``pressure_pa``: P_sat(T) = y_w·P.
    None when the gas holds no water.  Below 0 °C: over supercooled water.
    ValueError for invalid inputs or a dew point outside −50 … 373.85 °C."""
    y, p = water_mole_fraction_gas, pressure_pa
    if not (math.isfinite(y) and 0.0 <= y <= 1.0 and math.isfinite(p) and p > 0.0):
        raise ValueError(f"need a water mole fraction in [0, 1] and a pressure > 0 (got {y}, {p} Pa)")
    target = y * p
    if target <= 0.0:
        return None
    lo, hi = 223.15, 647.0                      # −50 °C … near the critical point
    if water_psat_pa(lo) > target:
        raise ValueError(f"water dew point below −50 °C (y_w = {y:.4g} at {p:.4g} Pa)")
    if water_psat_pa(hi) < target:
        raise ValueError(f"water dew point above {hi - 273.15:.2f} °C (y_w·P = {target:.4g} Pa)")
    for _ in range(100):                        # P_sat increases with T
        mid = 0.5 * (lo + hi)
        if water_psat_pa(mid) < target:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-7:
            break
    return 0.5 * (lo + hi) - 273.15
