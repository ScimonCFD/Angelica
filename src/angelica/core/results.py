from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class ComponentFlowResult:
    label: str
    mass_flow_kg_per_s: float
    volumetric_flow_m3_per_h: float
    temperature_in_c: float | None = None
    temperature_out_c: float | None = None


@dataclass(frozen=True)
class IterationMetrics:
    pressure_correction_abs_pa: float
    pressure_correction_mean_abs_pa: float
    pressure_correction_rel: float
    max_nodal_mass_imbalance_rel: float


@dataclass(frozen=True)
class SolveResult:
    case_name: str
    converged: bool
    node_pressures_pa: dict[int, float]
    component_flows: list[ComponentFlowResult]
    laminar_history: list[float]
    laminar_metrics: list[IterationMetrics]
    turbulent_history: list[float]
    turbulent_metrics: list[IterationMetrics]
    node_temperatures_c: dict[int, float] = field(default_factory=dict)
    temperature_history: list[float] = field(default_factory=list)
    density_history: list[float] = field(default_factory=list)
    outer_turbulent_final_metrics: tuple[IterationMetrics, ...] = field(default_factory=tuple)

    @property
    def link_mass_flows_kg_per_s(self) -> list[float]:
        return [component.mass_flow_kg_per_s for component in self.component_flows]

    @property
    def link_volumetric_flows_m3_per_h(self) -> list[float]:
        return [component.volumetric_flow_m3_per_h for component in self.component_flows]
