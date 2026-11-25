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


WELL_STATES = ("producing", "dead_column", "dead_back_pressure", "dead_minimum_rate", "unstable")


@dataclass(frozen=True)
class WellOperatingState:
    """Operating state of a ``Well`` (IPR + tubing) at the solution.

    ``state``:
      ``"producing"``           on the stable branch of its deliverability curve
      ``"dead_column"``         the static column outweighs the reservoir: it
                                can never flow without artificial lift
      ``"dead_back_pressure"``  wellhead above the curve's maximum
      ``"dead_minimum_rate"``   wellhead between the static pressure and the
                                maximum: producing would need at least the
                                minimum stable rate, which raises the network
                                above the maximum.  A dead well cannot start
                                against this pressure; unloading or gas lift
                                would be needed.
      ``"unstable"``            no steady state: dead, its wellhead is below
                                the static pressure, so it would start;
                                producing at its minimum stable rate raises
                                the wellhead above the maximum, so it would
                                die.  The real well produces intermittently.
                                The solve reports ``converged=False``.
    """

    component_id: str
    state: str
    wellhead_pressure_pa: float
    max_wellhead_pressure_pa: float      # maximum of the deliverability curve
    static_wellhead_pressure_pa: float   # curve at zero rate
    minimum_stable_rate_kg_per_s: float  # rate at the maximum
    mass_flow_kg_per_s: float            # reservoir production
    gas_lift_kg_per_s: float = 0.0       # lift gas injected (leaves at the wellhead)

    @property
    def message(self) -> str:
        mpa = 1e-6
        if self.state == "producing":
            return "producing"
        if self.state == "dead_column":
            return "dead: column heavier than P_r"
        if self.state == "dead_back_pressure":
            return f"dead: needs P_wh < {self.max_wellhead_pressure_pa * mpa:.2f} MPa"
        if self.state == "unstable":
            return (f"unstable, no steady state: dead it would start (P_wh "
                    f"{self.wellhead_pressure_pa * mpa:.2f} < static {self.static_wellhead_pressure_pa * mpa:.2f} MPa); "
                    f"producing its minimum stable rate ({self.minimum_stable_rate_kg_per_s:.2f} kg/s) "
                    f"it would exceed its maximum {self.max_wellhead_pressure_pa * mpa:.2f} MPa")
        return (f"dead: cannot sustain its minimum stable rate "
                f"({self.minimum_stable_rate_kg_per_s:.2f} kg/s) against "
                f"{self.wellhead_pressure_pa * mpa:.2f} MPa; to start it needs "
                f"P_wh < {self.static_wellhead_pressure_pa * mpa:.2f} MPa, unloading or gas lift")


@dataclass(frozen=True)
class SeparatorOutletResult:
    """One ``SeparatorOutlet`` link at the solution."""

    component_id: str
    phase: str
    mass_flow_kg_per_s: float
    valve_pressure_drop_pa: float    # P_separator − P_outlet node (level-control valve)


@dataclass(frozen=True)
class SeparatorResult:
    """Phase split of a separator node at the solution.

    Masses are of the phases in the vessel (gas, oil with its dissolved gas,
    free water).  Standard rates: the gas phase's gas, the oil phase's
    stock-tank liquid and the water, at 101 325 Pa and 15.56 °C.
    """

    node_id: int
    pressure_pa: float
    temperature_c: float
    feed_mass_flow_kg_per_s: float
    gas_mass_flow_kg_per_s: float
    oil_mass_flow_kg_per_s: float
    water_mass_flow_kg_per_s: float
    gas_rate_sm3_per_day: float
    oil_rate_sm3_per_day: float
    water_rate_m3_per_day: float
    outlets: tuple[SeparatorOutletResult, ...] = ()
    # Links through which fluid flows back into the separator from a node
    # with no declared fluid (e.g. a gas export above the separator
    # pressure).  That fluid is taken as the separator's own gas.
    backflow_links: tuple[str, ...] = ()

    @property
    def backflow_without_gas(self) -> bool:
        """Back-flow into a vessel with no gas phase: the inflow can only
        leave dissolved in the liquid, and how much would dissolve is not
        defined (no fluid declared upstream), so there is no steady state."""
        return bool(self.backflow_links) and self.gas_mass_flow_kg_per_s <= 1e-9 * max(
            self.feed_mass_flow_kg_per_s, 1e-30)

    @property
    def message(self) -> str:
        parts = []
        low = [o.component_id or o.phase for o in self.outlets if o.valve_pressure_drop_pa < 0.0]
        if low:
            parts.append("outlet pressure above the separator at " + ", ".join(low)
                         + ": the liquid needs a pump")
        if self.backflow_without_gas:
            parts.append("fluid flows back into the separator through " + ", ".join(self.backflow_links)
                         + " but the vessel has no gas phase: no steady state")
        elif self.backflow_links:
            parts.append("gas flows back into the separator through " + ", ".join(self.backflow_links)
                         + " (taken as the separator's own gas; declare a fluid at the boundary to set it)")
        return "; ".join(parts) if parts else "ok"


@dataclass(frozen=True)
class RegulatorResult:
    """A ``PressureRegulator`` at the solution.

    ``state``: "controlling" (the controlled node is at the set point),
    "open" (wide open, set point not reached) or "closed" (no flow).
    """

    component_id: str
    mode: str
    set_pressure_pa: float
    state: str
    inlet_pressure_pa: float
    outlet_pressure_pa: float
    mass_flow_kg_per_s: float

    @property
    def pressure_drop_pa(self) -> float:
        return self.inlet_pressure_pa - self.outlet_pressure_pa

    @property
    def message(self) -> str:
        held = "outlet" if self.mode == "downstream" else "inlet"
        if self.state == "controlling":
            return f"holding the {held} at the set point"
        if self.state == "open":
            return (f"wide open: the {held} cannot reach the set point "
                    + ("(upstream pressure too low)" if self.mode == "downstream"
                       else "(downstream pressure too high)"))
        return ("closed: the outlet is at or above the set point without it" if self.mode == "downstream"
                else "closed: the inlet is below the set point")


@dataclass(frozen=True)
class PatternSwitch:
    """A pipe whose Beggs-Brill flow pattern kept changing in the last outer
    iterations of a solve that did not converge."""

    component_id: str
    patterns: tuple[str, ...]


@dataclass(frozen=True)
class ConvergenceDiagnosis:
    """Where a solve that did not converge is stuck (``SolveResult.diagnosis``).

    ``failed``: the criteria not met ("hydraulics", "temperature",
    "density", "separator split", "wells", or a separator's stated reason).
    ``worst_node_id``: the interior node with the largest relative mass
    imbalance |in − out|/(in + out) at the end, the hydraulic convergence
    measure.  ``pattern_switches``: pipes whose flow pattern changed over the
    last outer iterations — a flow-pattern boundary the solution cannot settle
    on (Beggs-Brill's pattern boundaries are discontinuous, so a pipe can have
    no solution at an imposed flow, or two at one pressure drop).
    """

    failed: tuple[str, ...]
    worst_node_id: int | None = None
    worst_node_imbalance_rel: float = 0.0
    pattern_switches: tuple[PatternSwitch, ...] = ()

    @property
    def message(self) -> str:
        parts = ["not converged: " + ", ".join(self.failed)]
        if self.worst_node_id is not None and self.worst_node_imbalance_rel > 0.0:
            parts.append(f"largest mass imbalance at node {self.worst_node_id} "
                         f"({100.0 * self.worst_node_imbalance_rel:.2g} % of its throughput)")
        for sw in self.pattern_switches:
            parts.append(f"flow pattern of {sw.component_id} alternates between " + " and ".join(sw.patterns))
        return "; ".join(parts)


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
    well_states: tuple[WellOperatingState, ...] = ()
    separators: tuple[SeparatorResult, ...] = ()
    regulators: tuple[RegulatorResult, ...] = ()
    # None when converged (see ConvergenceDiagnosis)
    diagnosis: ConvergenceDiagnosis | None = None

    @property
    def link_mass_flows_kg_per_s(self) -> list[float]:
        return [component.mass_flow_kg_per_s for component in self.component_flows]

    @property
    def link_volumetric_flows_m3_per_h(self) -> list[float]:
        return [component.volumetric_flow_m3_per_h for component in self.component_flows]
