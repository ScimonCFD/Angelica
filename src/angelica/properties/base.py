from __future__ import annotations

from abc import ABC, abstractmethod


class FluidModel(ABC):
    @abstractmethod
    def density_for_link(self, link_state) -> float:
        raise NotImplementedError

    @abstractmethod
    def viscosity_for_link(self, link_state) -> float:
        raise NotImplementedError
