from __future__ import annotations

from abc import ABC, abstractmethod

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

    def density(self, pressure_pa: float, temperature_c: float) -> float:
        from thermo.eos import PR as _PR
        T = temperature_c + 273.15
        pr = _PR(
            Tc=self.critical_temperature_k,
            Pc=self.critical_pressure_pa,
            omega=self.acentric_factor,
            T=T,
            P=pressure_pa,
        )
        # Prefer the gas root; fall back to liquid root when only one root exists.
        Z = getattr(pr, "Z_g", None) or getattr(pr, "Z_l", None)
        if Z is None:
            raise RuntimeError(
                f"PR EOS: no valid Z root at T={temperature_c:.1f} °C, "
                f"P={pressure_pa / 1e6:.3f} MPa"
            )
        return pressure_pa * self.molecular_weight_kg_per_mol / (Z * _R * T)
