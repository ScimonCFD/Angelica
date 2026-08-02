from __future__ import annotations

from math import exp
from typing import Callable

from .eos import EquationOfState


def lee_gonzalez_eakin_viscosity(
    eos: EquationOfState,
    molecular_weight_kg_per_mol: float,
) -> Callable[[float, float], float]:
    """Natural gas viscosity via the Lee-Gonzalez-Eakin (1966) correlation.

    Returns a callable f(pressure_pa, temperature_c) → viscosity_pa_s
    suitable for use as the viscosity argument of CompressibleFluid.from_functions().

    The correlation evaluates local gas density from the supplied EOS, so
    it captures the pressure dependence of viscosity for real gases. At
    near-atmospheric conditions it converges to the temperature-only limit.

    Typical error < 2 % for natural gas in the range 38–171 °C, SG 0.52–0.80.

    Args:
        eos: Equation of state used to evaluate gas density at (P, T).
        molecular_weight_kg_per_mol: Molar mass of the gas (kg/mol).
    """
    Mg = molecular_weight_kg_per_mol * 1000.0  # kg/mol → g/mol

    def viscosity_pa_s(pressure_pa: float, temperature_c: float) -> float:
        T_R = (temperature_c + 273.15) * 1.8              # °C → Rankine
        rho_gcc = eos.density(pressure_pa, temperature_c) / 1000.0  # kg/m³ → g/cm³
        K = (9.4 + 0.02 * Mg) * T_R ** 1.5 / (209.0 + 19.0 * Mg + T_R)
        X = 3.5 + 986.0 / T_R + 0.01 * Mg
        Y = 2.4 - 0.2 * X
        mu_cP = K * exp(X * rho_gcc ** Y) * 1e-4
        return mu_cP * 1e-3  # cP → Pa·s

    return viscosity_pa_s
