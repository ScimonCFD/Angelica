from __future__ import annotations

from abc import ABC, abstractmethod


class FluidModel(ABC):
    @abstractmethod
    def density_for_link(self, link_state) -> float:
        raise NotImplementedError

    @abstractmethod
    def viscosity_for_link(self, link_state) -> float:
        raise NotImplementedError

    def specific_heat_for_link(self, _link_state) -> float:
        raise NotImplementedError(
            f"{type(self).__name__} does not provide specific heat capacity. "
            "Use ThermalFluid for non-isothermal simulations."
        )

    def thermal_conductivity_for_link(self, _link_state) -> float:
        raise NotImplementedError(
            f"{type(self).__name__} does not provide thermal conductivity. "
            "Use ThermalFluid for non-isothermal simulations."
        )
