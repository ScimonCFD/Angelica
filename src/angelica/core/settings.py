from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class SolverSettings:
    laminar_iterations: Optional[int] = None
    laminar_iterations_without_fittings: int = 1
    laminar_iterations_with_fittings: int = 7
    turbulent_iterations: int = 60
    pressure_relaxation: float = 1.0
    colebrook_friction_strategy: str = "transformed"
    friction_factor_method: str = "newton"
    friction_factor_max_iterations: int = 50
    velocity_loop_method: str = "fixed_point"
    velocity_loop_max_iterations: int = 50
    colebrook_residual_tolerance: float = 1e-4
    velocity_loop_tolerance: float = 1e-4
    pressure_correction_abs_tolerance_pa: float = 1e-3
    pressure_correction_rel_tolerance: float = 1e-8
    nodal_mass_imbalance_tolerance_kg_per_s: float = 1e-4
    nodal_mass_imbalance_rel_tolerance: float = 1e-3

    @property
    def friction_tolerance(self) -> float:
        return self.colebrook_residual_tolerance

    @property
    def correction_tolerance(self) -> float:
        return self.pressure_correction_abs_tolerance_pa
