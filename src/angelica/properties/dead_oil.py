from __future__ import annotations

import math


def dead_oil_specific_heat_j_per_kg_k(api_gravity: float, temperature_c: float) -> float:
    """Dead oil specific heat capacity via Watson-Nelson (API Technical Data Book).

    Correlation: cp [BTU/(lb·°F)] = (0.388 + 0.00045·T_F) / sqrt(SG)
    where T_F is temperature in Fahrenheit and SG is specific gravity at 60 °F.

    Applicable range: API 10–70, -18 °C to 260 °C.
    Returns specific heat in J/(kg·K).
    """
    specific_gravity = 141.5 / (api_gravity + 131.5)
    temperature_f = temperature_c * 9.0 / 5.0 + 32.0
    cp_btu = (0.388 + 0.00045 * temperature_f) / math.sqrt(specific_gravity)
    return cp_btu * 4186.8


def dead_oil_thermal_conductivity_w_per_m_k(api_gravity: float, temperature_c: float) -> float:
    """Dead oil thermal conductivity via Cragoe (1929).

    Correlation: k [BTU/(h·ft·°F)] = 0.0677·(1 - 0.0003·(T_F - 32)) / sqrt(SG)
    Converted to SI: k [W/(m·K)] = 0.1172·(1 - 0.00054·T_C) / sqrt(SG)

    Applicable range: API 10–70, up to ~200 °C.
    Returns thermal conductivity in W/(m·K).
    """
    specific_gravity = 141.5 / (api_gravity + 131.5)
    return 0.1172 * (1.0 - 0.00054 * temperature_c) / math.sqrt(specific_gravity)


def build_thermal_dead_oil(api_gravity: float):
    """Return a ThermalFluid with temperature-dependent dead oil properties.

    Density: Beggs & Robinson standard-condition density (constant w.r.t. T,
             consistent with incompressible assumption).
    Viscosity: Beggs & Robinson (1975) dead-oil correlation.
    Specific heat: Watson-Nelson (API Technical Data Book).
    Thermal conductivity: Cragoe (1929).

    No range checking is performed beyond the hard limit on API gravity in
    dead_oil_density_kg_per_m3. The caller is responsible for ensuring that
    operating temperatures and API values lie within the stated validity ranges
    of each underlying correlation.

    Args:
        api_gravity: API gravity of the dead oil (°API). Typical range 10–58.

    Returns:
        ThermalFluid instance ready for use with SteadyNonIsothermalIncompressibleSolver.
    """
    from .thermal_fluid import ThermalFluid

    rho = dead_oil_density_kg_per_m3(api_gravity)

    return ThermalFluid.from_functions(
        density_fn=lambda _t: rho,
        viscosity_fn=lambda t: dead_oil_viscosity_pa_s(api_gravity, t),
        specific_heat_fn=lambda t: dead_oil_specific_heat_j_per_kg_k(api_gravity, t),
        thermal_conductivity_fn=lambda t: dead_oil_thermal_conductivity_w_per_m_k(api_gravity, t),
    )


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

    Correlation stated range: 16 <= API <= 58, -6.7 °C <= T <= 146 °C.
    Hard lower bound: T > 0 °F (> -17.78 °C) — log10(0) is undefined below this point.
    Returns viscosity in Pa·s.
    """
    temperature_f = temperature_c * 9.0 / 5.0 + 32.0
    if temperature_f <= 0.0:
        raise ValueError(f"Temperature must be above 0 °F (-17.78 °C); got {temperature_f:.1f} °F ({temperature_c:.1f} °C).")
    z = 3.0324 - 0.02023 * api_gravity - 1.163 * math.log10(temperature_f)
    viscosity_cp = 10.0 ** (10.0 ** z) - 1.0
    return max(viscosity_cp / 1000.0, 1e-9)
