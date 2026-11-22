"""
Gas-liquid surface tension for black-oil fluids.

Beggs-Brill needs the liquid surface tension only through the liquid velocity
number N_LV = v_sL·(ρ_L/(g·σ))^¼, so the correlations below are the usual
black-oil ones (Brill & Mukherjee, *Multiphase Flow in Wells*, SPE Monograph
17, 1999, Appendix B; Beggs, *Production Optimization Using Nodal Analysis*,
1991):

Gas–oil, Baker & Swerdloff (1956)
    dead oil  σ_68  = 39.0 − 0.2571·API   (dyn/cm, at 68 °F)
              σ_100 = 37.5 − 0.2571·API   (dyn/cm, at 100 °F)
              linear in T between 68 and 100 °F, constant outside
    live oil  σ_o = σ_od·(1 − 0.024·p^0.45),  p in psia

Gas–water (Hough et al. curves as fitted in the references above)
    σ_w(74 °F)  = 75 − 1.108·p^0.349
    σ_w(280 °F) = 53 − 0.1048·p^0.637,  p in psia, dyn/cm
    linear in T between 74 and 280 °F, constant outside

Liquid surface tension = no-slip volume average of oil and water.
All values are floored at 1 dyn/cm (the live-oil factor reaches zero near
4 000 psia, where the correlation is no longer meaningful).
"""
from __future__ import annotations

_PSI_PER_PA = 1.0 / 6894.757
_DYN_PER_CM_TO_N_PER_M = 1.0e-3
_FLOOR_DYN_PER_CM = 1.0


def _fahrenheit(t_c: float) -> float:
    return t_c * 9.0 / 5.0 + 32.0


def dead_oil_surface_tension(api_gravity: float, temperature_c: float) -> float:
    """Baker-Swerdloff dead-oil surface tension (N/m)."""
    s68 = 39.0 - 0.2571 * api_gravity
    s100 = 37.5 - 0.2571 * api_gravity
    t_f = _fahrenheit(temperature_c)
    if t_f <= 68.0:
        s = s68
    elif t_f >= 100.0:
        s = s100
    else:
        s = s68 - (t_f - 68.0) * (s68 - s100) / 32.0
    return max(s, _FLOOR_DYN_PER_CM) * _DYN_PER_CM_TO_N_PER_M


def oil_surface_tension(api_gravity: float, temperature_c: float, pressure_pa: float) -> float:
    """Gas–oil surface tension of live oil (N/m), Baker-Swerdloff."""
    p_psia = max(pressure_pa, 0.0) * _PSI_PER_PA
    factor = 1.0 - 0.024 * p_psia ** 0.45
    s = dead_oil_surface_tension(api_gravity, temperature_c) / _DYN_PER_CM_TO_N_PER_M * factor
    return max(s, _FLOOR_DYN_PER_CM) * _DYN_PER_CM_TO_N_PER_M


def water_surface_tension(temperature_c: float, pressure_pa: float) -> float:
    """Gas–water surface tension (N/m)."""
    p_psia = max(pressure_pa, 0.0) * _PSI_PER_PA
    s74 = 75.0 - 1.108 * p_psia ** 0.349
    s280 = 53.0 - 0.1048 * p_psia ** 0.637
    t_f = _fahrenheit(temperature_c)
    if t_f <= 74.0:
        s = s74
    elif t_f >= 280.0:
        s = s280
    else:
        s = s74 - (t_f - 74.0) * (s74 - s280) / 206.0
    return max(s, _FLOOR_DYN_PER_CM) * _DYN_PER_CM_TO_N_PER_M


def liquid_surface_tension(
    sigma_oil: float, sigma_water: float, oil_volume: float, water_volume: float
) -> float:
    """No-slip volume average of the oil and water surface tensions (N/m)."""
    total = oil_volume + water_volume
    if total <= 0.0:
        return sigma_oil
    return (oil_volume * sigma_oil + water_volume * sigma_water) / total
