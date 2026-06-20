from dataclasses import dataclass

from .base import FluidModel


@dataclass(frozen=True)
class SingleComponentFluid(FluidModel):
    density_kg_per_m3: float
    viscosity_pa_s: float

    def __post_init__(self) -> None:
        if self.density_kg_per_m3 <= 0.0:
            raise ValueError(f"density_kg_per_m3 must be positive (got {self.density_kg_per_m3})")
        if self.viscosity_pa_s <= 0.0:
            raise ValueError(f"viscosity_pa_s must be positive (got {self.viscosity_pa_s})")

    def density_for_link(self, link_state) -> float:
        return self.density_kg_per_m3

    def viscosity_for_link(self, link_state) -> float:
        return self.viscosity_pa_s
