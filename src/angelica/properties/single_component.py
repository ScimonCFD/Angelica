from dataclasses import dataclass

from .base import FluidModel


@dataclass(frozen=True)
class SingleComponentFluid(FluidModel):
    density_kg_per_m3: float
    viscosity_pa_s: float

    def density_for_link(self, link_state) -> float:
        return self.density_kg_per_m3

    def viscosity_for_link(self, link_state) -> float:
        return self.viscosity_pa_s
