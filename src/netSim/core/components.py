from __future__ import annotations

from dataclasses import dataclass
import math


@dataclass(frozen=True)
class PressureChanger:
    start_node: int
    end_node: int
    diameter_m: float
    component_id: str = ""

    @property
    def area_m2(self) -> float:
        return math.pi * self.diameter_m**2 / 4.0


@dataclass(frozen=True)
class Pipe(PressureChanger):
    length_m: float = 0.0
    absolute_roughness_m: float = 0.0
    hazen_williams_c: float = 130.0
    height_change_m: float = 0.0


@dataclass(frozen=True)
class Fitting(PressureChanger):
    loss_coefficient: float = 0.0
