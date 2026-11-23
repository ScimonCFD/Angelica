from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True)
class MultiphaseFlowResult:
    """Two-phase flow state of one pipe (black-oil and compositional runs).

    ``flow_pattern`` is the Beggs-Brill pattern when the run used Beggs-Brill,
    or ``"no-slip"`` for the homogeneous model (then the holdup equals the
    no-slip liquid fraction).  Gradients are along the flow direction, positive
    when they oppose the flow; the elevation gradient uses the in-situ holdup.
    The erosional velocity is API RP 14E's with C = 100 (continuous service)
    on the no-slip mixture density.
    """

    flow_pattern: str
    no_slip_liquid_fraction: float          # C_L = q_L / (q_L + q_g)
    liquid_holdup: float                    # H_L, in-situ liquid volume fraction
    superficial_liquid_velocity_m_per_s: float
    superficial_gas_velocity_m_per_s: float
    mixture_velocity_m_per_s: float
    liquid_volume_m3: float                 # H_L · A · L
    elevation_gradient_pa_per_m: float
    friction_gradient_pa_per_m: float
    erosional_velocity_m_per_s: float = float("nan")   # API RP 14E, C = 100
    erosional_velocity_ratio: float = float("nan")      # v_m / V_e (> 1: above the limit)


@dataclass(frozen=True)
class ComponentFlowResult:
    label: str
    mass_flow_kg_per_s: float
    volumetric_flow_m3_per_h: float
    temperature_in_c: float | None = None
    temperature_out_c: float | None = None
    zs: tuple[float, ...] = ()        # mole fractions — non-empty for compositional runs
    vapor_fraction: float | None = None  # 0=liquid, 1=gas, (0,1)=two-phase; None for non-EOS modes
    free_water_fraction: float = 0.0   # compositional: mole fraction of the feed that is liquid free water;
                                       # black-oil: in-situ water volume fraction
    liquid_fraction: float | None = None  # mole fraction of feed that is HC liquid
    gas_phase_zs: tuple[float, ...] = ()    # gas phase mole fractions per HC component
    liquid_phase_zs: tuple[float, ...] = () # liquid HC phase mole fractions per HC component
    multiphase: MultiphaseFlowResult | None = None  # pipes in black-oil / compositional runs
    # Hydrate check (compositional pipes with hydrate formers), at the pipe end
    # where the fluid is coldest relative to the formation curve:
    hydrate_temperature_c: float | None = None   # formation temperature at that end's pressure
    hydrate_subcooling_k: float | None = None    # T_hydrate − T_fluid (> 0: inside the hydrate region)
    hydrate_structure: str = ""


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
        shaft_work_kw    > 0 : work ADDED by pumps / compressors
        joule_thomson_kw     : change of ṁ·cp·T at constant enthalpy across
                               valves and along pipes (< 0 when a real gas
                               cools on expansion); part of the balance only
                               with enthalpy_basis = "sensible"
        potential_energy_kw  : −Σ ṁ·g·Δz along pipes (< 0 when lifting)
        energy_error_kw ≈ 0  : residual of
                               (Ė_in + Q_src + W_shaft + Ė_JT + Ė_pot − Q_wall − Ė_out)
    """

    enthalpy_in_kw: float       # Σ ṁ × cp × T at inlet boundary nodes
    enthalpy_out_kw: float      # Σ ṁ × cp × T at outlet boundary nodes
    heat_sources_kw: float      # Σ power_w for all HeatSource links
    heat_wall_loss_kw: float    # Σ U·π·D·L·(T_pipe − T_amb) for all pipes + heat wells give to the formation
    shaft_work_kw: float = 0.0  # Σ pump / compressor work on the fluid
    joule_thomson_kw: float = 0.0  # Σ ṁ·Δ(cp·T) at constant h (valves and pipes)
    potential_energy_kw: float = 0.0  # −Σ ṁ·g·Δz along pipes and well tubing (enthalpy-based fluids)
    # "sensible": fluxes are ṁ·cp·T (liquids, black-oil); "real": fluxes are the
    # fluid's enthalpy ṁ·h(P, T) above the ideal gas at 0 °C (compressible gas,
    # compositional EOS).  With real enthalpy the Joule-Thomson effect is part
    # of h, so joule_thomson_kw is reported for information only.
    enthalpy_basis: str = "sensible"

    @property
    def energy_error_kw(self) -> float:
        """Residual of the global first-law balance (should be ≈ 0)."""
        jt = self.joule_thomson_kw if self.enthalpy_basis == "sensible" else 0.0
        return (self.enthalpy_in_kw + self.heat_sources_kw + self.shaft_work_kw
                + jt + self.potential_energy_kw
                - self.heat_wall_loss_kw - self.enthalpy_out_kw)

    @property
    def energy_error_pct(self) -> float:
        """Residual as a percentage of the largest flux.  With real enthalpy the
        boundary fluxes depend on the reference state (they can be small next
        to the heat exchanged or the Joule-Thomson effect), so every term of
        the balance is considered.

        The two bases are not comparable: for the same residual in kW, the
        "real" percentage (largest of all terms) is usually smaller than a
        "sensible" one (largest boundary flux).  Judge accuracy by
        ``energy_error_kw`` against the heat of a small temperature change of
        the flow, e.g. ṁ·cp·0.01 K."""
        ref = max(abs(self.enthalpy_in_kw), abs(self.enthalpy_out_kw), 1e-30)
        if self.enthalpy_basis == "real":
            ref = max(ref, abs(self.heat_sources_kw), abs(self.heat_wall_loss_kw),
                      abs(self.shaft_work_kw), abs(self.joule_thomson_kw), abs(self.potential_energy_kw))
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
    density_rel_change_history: list[float] = field(default_factory=list)

    @property
    def density_history(self) -> "list[float]":
        """Backward-compatible alias for density_rel_change_history."""
        return self.density_rel_change_history
    outer_turbulent_final_metrics: tuple[IterationMetrics, ...] = field(default_factory=tuple)
    outer_iteration_boundaries: tuple[int, ...] = field(default_factory=tuple)
    global_balance: GlobalBalance | None = None
    global_energy_balance: GlobalEnergyBalance | None = None
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
