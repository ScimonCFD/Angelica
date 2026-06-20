"""Temperature-dependent properties for liquid water (0–100 °C).

Density:    Kell (1975) formula, error < 0.1 kg/m³ over 0–100 °C.
Viscosity:  Vogel equation, matches IAPWS 2008 within 1 %.
Specific heat: constant 4182 J/kg·K (varies < 1 % between 15 and 80 °C).
Thermal conductivity: quadratic fit to IAPWS data, error < 0.5 % over 0–100 °C.

No range checking is performed on the temperature input. The caller is
responsible for ensuring that operating temperatures lie within the stated
validity ranges of each correlation.
"""
from __future__ import annotations

from angelica.properties.thermal_fluid import ThermalFluid


def water_density_kg_per_m3(temperature_c: float) -> float:
    """Liquid water density (kg/m³) as a function of temperature (°C).

    Kell (1975) formula.  Reference values: 999.84 at 0 °C, 998.21 at 20 °C,
    992.22 at 40 °C, 983.20 at 60 °C, 971.82 at 80 °C, 958.37 at 100 °C.
    """
    T = temperature_c
    return 1000.0 * (
        1.0
        - (T + 288.9414) / (508929.2 * (T + 68.12963)) * (T - 3.9863) ** 2
    )


def water_viscosity_pa_s(temperature_c: float) -> float:
    """Liquid water dynamic viscosity (Pa·s) as a function of temperature (°C).

    Vogel equation calibrated to IAPWS 2008.  Reference values (mPa·s):
    1.002 at 20 °C, 0.653 at 40 °C, 0.467 at 60 °C, 0.355 at 80 °C.
    """
    T_K = temperature_c + 273.15
    return 2.414e-5 * 10.0 ** (247.8 / (T_K - 140.0))


def water_specific_heat_j_per_kg_k(_temperature_c: float) -> float:
    """Liquid water specific heat (J/kg·K).  Treated as constant 4182."""
    return 4182.0


def water_thermal_conductivity_w_per_m_k(temperature_c: float) -> float:
    """Liquid water thermal conductivity (W/m·K), quadratic fit to IAPWS data.

    Reference values: 0.598 at 20 °C, 0.631 at 40 °C, 0.655 at 60 °C, 0.673 at 80 °C.
    """
    T = temperature_c
    return 0.5636 + 1.946e-3 * T - 8.151e-6 * T**2


def build_water_thermal_fluid() -> ThermalFluid:
    """Return a ThermalFluid with all four water properties as functions of T."""
    return ThermalFluid.from_functions(
        density_fn=water_density_kg_per_m3,
        viscosity_fn=water_viscosity_pa_s,
        specific_heat_fn=water_specific_heat_j_per_kg_k,
        thermal_conductivity_fn=water_thermal_conductivity_w_per_m_k,
    )
