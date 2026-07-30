from __future__ import annotations

from abc import ABC, abstractmethod

_GAS_CONSTANT_J_PER_MOL_K = 8.314


class EquationOfState(ABC):
    @abstractmethod
    def density(self, pressure_pa: float, temperature_c: float) -> float:
        raise NotImplementedError


class IdealGasEOS(EquationOfState):
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
            / (_GAS_CONSTANT_J_PER_MOL_K * (temperature_c + 273.15))
        )
