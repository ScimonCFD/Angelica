"""Free-water (immiscible-water) utilities for compositional flash.

When water is present in the component list the compositional solver separates
it from the hydrocarbon system, performs the VL flash on the HC-only
normalised composition, and then determines how much water is liquid or vapour
using the water saturation pressure (Wagner / IAPWS-IF97 correlation).

This module provides only the pure physics:
  - water detection in a component list
  - water saturation pressure (Wagner equation, IAPWS-IF97 accuracy)
  - liquid water density (Kell 1975 correlation)
  - liquid water dynamic viscosity (empirical fit)

The actual flash coupling lives in ``compositional_fluid.py``.
"""

from __future__ import annotations

import math

# ── water identifiers ────────────────────────────────────────────────────────

_WATER_NAMES_LOWER: frozenset[str] = frozenset(
    {"water", "h2o", "agua", "7732-18-5", "oxidane", "dihydrogen oxide", "dihydrogen monoxide"}
)

WATER_MW: float = 18.01528  # g/mol
_TC: float = 647.096  # K  (IAPWS-IF97 critical temperature)
_PC: float = 22_064_000.0  # Pa (IAPWS-IF97 critical pressure)


def find_water_index(component_names: tuple[str, ...]) -> int | None:
    """Return the index of water in the component list, or None."""
    for i, name in enumerate(component_names):
        if name.lower().strip() in _WATER_NAMES_LOWER:
            return i
    return None


# ── water saturation pressure (Wagner / IAPWS-IF97) ─────────────────────────

def water_psat_pa(T_K: float) -> float:
    """Water saturation pressure (Pa) via the Wagner equation.

    Coefficients from IAPWS-IF97.  Valid 273–647 K (0 °C to critical point).
    Error < 0.005 % throughout the range.
    """
    if T_K >= _TC:
        return _PC
    if T_K < 200.0:
        T_K = 200.0  # clamp — below 200 K water is solid in real conditions
    tau = 1.0 - T_K / _TC
    ln_pr = (_TC / T_K) * (
        -7.85951783 * tau
        + 1.84408259 * tau ** 1.5
        - 11.78664970 * tau ** 3
        + 22.68074110 * tau ** 3.5
        - 15.96187190 * tau ** 4
        + 1.80122502 * tau ** 7.5
    )
    return _PC * math.exp(ln_pr)


# ── liquid water properties ──────────────────────────────────────────────────

def water_liquid_density_kg_m3(T_K: float) -> float:
    """Liquid water density (kg/m³) via the Kell (1975) correlation.

    Valid 0–150 °C.  Clamped above 150 °C where the correlation extrapolates.
    """
    T_C = min(T_K - 273.15, 150.0)
    T_C = max(T_C, 0.0)
    num = (
        999.83952
        + 16.945176 * T_C
        - 7.987040e-3 * T_C ** 2
        - 46.170461e-6 * T_C ** 3
        + 105.56302e-9 * T_C ** 4
        - 280.54253e-12 * T_C ** 5
    )
    den = 1.0 + 16.879850e-3 * T_C
    return max(num / den, 700.0)


def water_liquid_viscosity_pa_s(T_K: float) -> float:
    """Liquid water dynamic viscosity (Pa·s).

    Vogel–Fulcher–Tammann fit.  Valid ~0–100 °C; extrapolates above that.
    """
    T_K = max(T_K, 273.15)
    return 2.414e-5 * 10.0 ** (247.8 / (T_K - 140.0))


# ── free-water mole-fraction split ──────────────────────────────────────────

def free_water_split(
    z_water: float,
    VF_hc: float,
    sum_hc: float,
    T_K: float,
    P_Pa: float,
) -> tuple[float, float]:
    """Return (n_water_liquid, n_water_vapour) mole fractions of the total feed.

    Uses the immiscible-water model: water forms a separate liquid phase
    when the partial pressure of water vapour in the gas stream would exceed
    the saturation pressure.

    The gas stream carries VF_hc * sum_hc moles of HC vapour per mole of
    feed.  Water vapour is in equilibrium with the gas phase at Raoult's law
    (P_water_vapour = y_w * P = Psat(T)):

        n_wv / (n_wv + V_HC) = Psat / P
        n_wv = V_HC * Psat / (P − Psat)

    Liquid water = z_water − n_wv (clamped to ≥ 0).

    Args:
        z_water:  mole fraction of water in the feed.
        VF_hc:    vapour fraction of the HC-only flash (0–1).
        sum_hc:   total HC mole fraction in the feed (= 1 − z_water).
        T_K:      temperature (K).
        P_Pa:     pressure (Pa).

    Returns:
        (n_water_liquid, n_water_vapour) — both ≥ 0, sum = z_water.
    """
    P_sat = water_psat_pa(T_K)
    V_HC = VF_hc * sum_hc  # HC vapour moles per mole of feed

    if P_Pa <= P_sat or V_HC < 1e-15:
        # Pressure below water saturation → all water stays as vapour
        return 0.0, z_water

    # Maximum water that can stay in the vapour phase
    max_n_wv = V_HC * P_sat / (P_Pa - P_sat)
    n_wv = min(z_water, max_n_wv)
    n_wl = z_water - n_wv
    return n_wl, n_wv
