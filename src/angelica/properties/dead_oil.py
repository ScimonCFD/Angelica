from __future__ import annotations

import math


def dead_oil_density_kg_per_m3(api_gravity: float) -> float:
    """Dead oil density at standard conditions (15.56 °C / 60 °F).

    Uses the API gravity definition: SG = 141.5 / (API + 131.5).
    Valid range: API > -131.5 (physical constraint).
    """
    if api_gravity <= -131.5:
        raise ValueError(f"API gravity must be greater than -131.5 (got {api_gravity}).")
    specific_gravity = 141.5 / (api_gravity + 131.5)
    return specific_gravity * 999.064


def dead_oil_viscosity_pa_s(api_gravity: float, temperature_c: float) -> float:
    """Dead oil dynamic viscosity via Beggs & Robinson (1975).

    Correlation: log10(log10(mu_od + 1)) = (3.0324 - 0.02023*API) - 1.163*log10(T_F)
    where T_F is temperature in Fahrenheit and mu_od is in centipoise.

    Stated range: 16 <= API <= 58, -6.7 °C <= T <= 146 °C.
    Returns viscosity in Pa·s.
    """
    temperature_f = temperature_c * 9.0 / 5.0 + 32.0
    if temperature_f <= 0.0:
        raise ValueError(f"Temperature must be above 0 °F; got {temperature_f:.1f} °F ({temperature_c:.1f} °C).")
    z = 3.0324 - 0.02023 * api_gravity - 1.163 * math.log10(temperature_f)
    viscosity_cp = 10.0 ** (10.0 ** z) - 1.0
    return max(viscosity_cp / 1000.0, 1e-9)
