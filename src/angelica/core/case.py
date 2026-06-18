from __future__ import annotations

from dataclasses import dataclass, field

from .components import PressureChanger
from angelica.properties.base import FluidModel


@dataclass(frozen=True)
class PressureBoundary:
    node_id: int
    pressure_pa: float


@dataclass(frozen=True)
class FlowBoundary:
    node_id: int
    mass_flow_kg_per_s: float


@dataclass(frozen=True)
class ThermalBoundary:
    node_id: int
    temperature_c: float = 0.0
    bc_type: str = "fixed_temperature"  # "fixed_temperature" | "zero_gradient" | "fixed_gradient"
    gradient_dc_per_m: float = 0.0


@dataclass(frozen=True)
class NetworkCase:
    name: str
    fluid_model: FluidModel
    pressure_inlets: tuple[PressureBoundary, ...]
    pressure_outlets: tuple[PressureBoundary, ...]
    components: tuple[PressureChanger, ...]
    flow_inlets: tuple[FlowBoundary, ...] = field(default_factory=tuple)
    flow_outlets: tuple[FlowBoundary, ...] = field(default_factory=tuple)
    node_ids: tuple[int, ...] = field(default_factory=tuple)
    initial_node_pressures_pa: dict[int, float] = field(default_factory=dict)
    thermal_inlets: tuple[ThermalBoundary, ...] = field(default_factory=tuple)

    def all_node_ids(self) -> tuple[int, ...]:
        if self.node_ids:
            return self.node_ids

        nodes: set[int] = set()
        for boundary in (
            self.pressure_inlets
            + self.pressure_outlets
            + self.flow_inlets
            + self.flow_outlets
        ):
            nodes.add(boundary.node_id)
        for component in self.components:
            nodes.add(component.start_node)
            nodes.add(component.end_node)
        return tuple(sorted(nodes))
