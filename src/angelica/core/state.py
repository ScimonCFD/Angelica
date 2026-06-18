from __future__ import annotations

from dataclasses import dataclass

from .components import Fitting, HeatSource, Pipe, PressureChanger, Pump


@dataclass
class NodeState:
    node_id: int
    pressure_pa: float | None = None
    is_boundary: bool = False
    is_inlet: bool = False
    is_pressure_boundary: bool = False
    prescribed_mass_flow_kg_per_s: float | None = None
    temperature_c: float | None = None
    is_thermal_inlet: bool = False
    thermal_gradient_dc_per_m: float | None = None


@dataclass
class PressureChangerState:
    component: PressureChanger
    start_node: NodeState
    end_node: NodeState
    velocity_m_per_s: float = 0.0
    mass_flow_kg_per_s: float = 0.0
    coupling: float = 0.0

    @property
    def area_m2(self) -> float:
        return self.component.area_m2

    @property
    def diameter_m(self) -> float:
        return self.component.diameter_m


@dataclass
class PipeState(PressureChangerState):
    component: Pipe
    friction_factor: float | None = None
    reynolds: float = 0.0
    temperature_c: float = 20.0


@dataclass
class FittingState(PressureChangerState):
    component: Fitting
    reynolds: float = 0.0


@dataclass
class PumpState(PressureChangerState):
    component: Pump
    reynolds: float = 0.0


@dataclass
class HeatSourceState(PressureChangerState):
    component: HeatSource
    reynolds: float = 0.0
    temperature_c: float = 20.0


@dataclass
class NetworkState:
    nodes: dict[int, NodeState]
    components: list[PressureChangerState]
