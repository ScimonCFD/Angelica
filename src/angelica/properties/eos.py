from __future__ import annotations

from abc import ABC, abstractmethod

import numpy as np

_R = 8.314  # J/(mol·K)


class EquationOfState(ABC):
    @abstractmethod
    def density(self, pressure_pa: float, temperature_c: float) -> float:
        raise NotImplementedError


class IdealGasEOS(EquationOfState):
    """Ideal gas: ρ = PM / RT."""

    def __init__(self, molecular_weight_kg_per_mol: float) -> None:
        if molecular_weight_kg_per_mol <= 0.0:
            raise ValueError(
                f"molecular_weight_kg_per_mol must be positive (got {molecular_weight_kg_per_mol})"
            )
        self.molecular_weight_kg_per_mol = molecular_weight_kg_per_mol

    def density(self, pressure_pa: float, temperature_c: float) -> float:
        return (
            pressure_pa
            * self.molecular_weight_kg_per_mol
            / (_R * (temperature_c + 273.15))
        )


class PengRobinsonEOS(EquationOfState):
    """Peng–Robinson cubic EOS for single-component gas.

    Solves the cubic compressibility-factor equation and returns the
    vapour-phase (largest real) root, giving density ρ = PM/(ZRT).
    Suitable for gas pipeline conditions away from the critical point
    and the two-phase region.

    Args:
        molecular_weight_kg_per_mol: Molar mass M (kg/mol).
        critical_temperature_k: Critical temperature Tc (K).
        critical_pressure_pa: Critical pressure Pc (Pa).
        acentric_factor: Pitzer acentric factor ω (dimensionless).
    """

    def __init__(
        self,
        molecular_weight_kg_per_mol: float,
        critical_temperature_k: float,
        critical_pressure_pa: float,
        acentric_factor: float,
    ) -> None:
        if molecular_weight_kg_per_mol <= 0.0:
            raise ValueError(f"molecular_weight_kg_per_mol must be positive (got {molecular_weight_kg_per_mol})")
        if critical_temperature_k <= 0.0:
            raise ValueError(f"critical_temperature_k must be positive (got {critical_temperature_k})")
        if critical_pressure_pa <= 0.0:
            raise ValueError(f"critical_pressure_pa must be positive (got {critical_pressure_pa})")

        self.molecular_weight_kg_per_mol = molecular_weight_kg_per_mol
        self.critical_temperature_k = critical_temperature_k
        self.critical_pressure_pa = critical_pressure_pa
        self.acentric_factor = acentric_factor

        kappa = 0.37464 + 1.54226 * acentric_factor - 0.26992 * acentric_factor ** 2
        self._kappa = kappa
        self._b = 0.07780 * _R * critical_temperature_k / critical_pressure_pa
        self._a_c = 0.45724 * _R ** 2 * critical_temperature_k ** 2 / critical_pressure_pa

    def _a(self, temperature_k: float) -> float:
        alpha = (1.0 + self._kappa * (1.0 - (temperature_k / self.critical_temperature_k) ** 0.5)) ** 2
        return self._a_c * alpha

    def density(self, pressure_pa: float, temperature_c: float) -> float:
        T = temperature_c + 273.15
        a = self._a(T)
        b = self._b
        A = a * pressure_pa / (_R * T) ** 2
        B = b * pressure_pa / (_R * T)
        # Z³ - (1-B)Z² + (A - 3B² - 2B)Z - (AB - B² - B³) = 0
        coeffs = [1.0, -(1.0 - B), A - 3.0 * B ** 2 - 2.0 * B, -(A * B - B ** 2 - B ** 3)]
        roots = np.roots(coeffs)
        # Vapour root: largest real root greater than B (the co-volume lower bound)
        real_roots = [
            r.real for r in roots
            if abs(r.imag) / max(abs(r.real), 1.0) < 1e-8 and r.real > B
        ]
        Z = max(real_roots) if real_roots else max(r.real for r in roots)
        return pressure_pa * self.molecular_weight_kg_per_mol / (Z * _R * T)
