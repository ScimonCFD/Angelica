from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional


@dataclass(frozen=True)
class ComponentFlowResult:
    label: str
    mass_flow_kg_per_s: float
    volumetric_flow_m3_per_h: float
    temperature_in_c: float | None = None
    temperature_out_c: float | None = None
    zs: tuple[float, ...] = ()        # mole fractions — non-empty for compositional runs
    vapor_fraction: float | None = None  # 0=liquid, 1=gas, (0,1)=two-phase; None for non-EOS modes


@dataclass(frozen=True)
class IterationMetrics:
    pressure_correction_abs_pa: float
    pressure_correction_mean_abs_pa: float
    pressure_correction_rel: float
    max_nodal_mass_imbalance_rel: float
    global_mass_imbalance_kg_per_s: float = 0.0
    global_mass_imbalance_rel: float = 0.0


@dataclass(frozen=True)
class GlobalBalance:
    mass_inlet_kg_per_s: float
    mass_outlet_kg_per_s: float

    @property
    def mass_error_pct(self) -> float:
        ref = max(self.mass_inlet_kg_per_s, self.mass_outlet_kg_per_s, 1e-30)
        return 100.0 * abs(self.mass_inlet_kg_per_s - self.mass_outlet_kg_per_s) / ref


@dataclass(frozen=True)
class GlobalEnergyBalance:
    """Steady-state global energy balance for thermal solvers.

    All quantities are in kW, relative to a 0 °C enthalpy reference.

    Sign convention:
        heat_sources_kw  > 0 : heat ADDED to the fluid (heaters, positive power_w)
        heat_wall_loss_kw > 0 : heat LEAVING the fluid to the environment
        energy_error_kw ≈ 0  : residual of (Ė_in + Q_src − Q_wall − Ė_out)
    """

    enthalpy_in_kw: float       # Σ ṁ × cp × T at inlet boundary nodes
    enthalpy_out_kw: float      # Σ ṁ × cp × T at outlet boundary nodes
    heat_sources_kw: float      # Σ power_w for all HeatSource links
    heat_wall_loss_kw: float    # Σ U·π·D·L·(T_pipe − T_amb) for all pipes

    @property
    def energy_error_kw(self) -> float:
        """Residual of the global first-law balance (should be ≈ 0)."""
        return self.enthalpy_in_kw + self.heat_sources_kw - self.heat_wall_loss_kw - self.enthalpy_out_kw

    @property
    def energy_error_pct(self) -> float:
        ref = max(abs(self.enthalpy_in_kw), abs(self.enthalpy_out_kw), 1e-30)
        return 100.0 * abs(self.energy_error_kw) / ref


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
    outer_iteration_boundaries: tuple[int, ...] = field(default_factory=tuple)
    global_balance: Optional[GlobalBalance] = None
    global_energy_balance: Optional[GlobalEnergyBalance] = None
    mass_balance_history: list[float] = field(default_factory=list)
    energy_balance_history: list[float] = field(default_factory=list)
    component_names: tuple[str, ...] = ()
    component_mws: tuple[float, ...] = ()
    node_compositions: dict[int, tuple[float, ...]] = field(default_factory=dict)

    @property
    def link_mass_flows_kg_per_s(self) -> list[float]:
        return [component.mass_flow_kg_per_s for component in self.component_flows]

    @property
    def link_volumetric_flows_m3_per_h(self) -> list[float]:
        return [component.volumetric_flow_m3_per_h for component in self.component_flows]
