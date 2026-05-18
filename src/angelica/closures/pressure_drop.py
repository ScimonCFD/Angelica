from __future__ import annotations

from abc import ABC, abstractmethod


class PressureDropCorrelation(ABC):
    @abstractmethod
    def calculate_velocity(
        self,
        link_state,
        delta_p: float,
        density: float,
        viscosity: float,
        tolerance: float | None = None,
    ) -> float:
        raise NotImplementedError

    @abstractmethod
    def calculate_coupling(self, link_state, density: float, viscosity: float) -> float:
        raise NotImplementedError
