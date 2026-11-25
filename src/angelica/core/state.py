from __future__ import annotations

from dataclasses import dataclass, field

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
    # Internal mass source (kg/s), e.g. lift gas leaving a gas-lifted well at
    # its wellhead node; ``source_streams`` describes it for composition
    # propagation (solver-specific descriptors).
    source_mass_flow_kg_per_s: float = 0.0
    source_streams: tuple = ()
    # Pressure held at a regulator's set point (an interior node whose
    # pressure is fixed while its ``PressureRegulator`` is controlling).
    regulated: bool = False


@dataclass
class PressureChangerState:
    component: PressureChanger
    start_node: NodeState
    end_node: NodeState
    velocity_m_per_s: float = 0.0
    mass_flow_kg_per_s: float = 0.0
    coupling: float = 0.0
    zs: tuple[float, ...] = ()        # mole fractions carried by the link (compositional runs)

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
    temperature_c: float | None = None
    exit_temperature_c: float | None = None   # fluid temperature leaving the pipe (flow direction)
    # Enthalpy-balance sources set by the energy solve (W, + heats the fluid):
    jt_source_w: float = 0.0        # −ṁ·[h(T̄, P_out) − h(T̄, P_in)]: Joule-Thomson
    gravity_source_w: float = 0.0   # −ṁ·g·Δz: potential energy (flow direction)


@dataclass
class FittingState(PressureChangerState):
    component: Fitting
    reynolds: float = 0.0
    # ``Well`` links only.  Between the static wellhead pressure and the
    # maximum of the deliverability curve a well is bistable (dead, or
    # producing on the stable branch); this remembers which state it is in.
    well_dead: bool = False
    well_latch_enabled: bool = True   # False: always produce when it can (P_wh < P_max)
    well_starts: int = 0              # times the solver restarted it (capped)
    # ``SeparatorOutlet`` links only: mass fraction of the separator feed that
    # leaves through this link, and the other links at the separator node
    # (their flows into the node are the feed).
    separator_fraction: float = 0.0
    separator_feed_links: tuple = ()
    # ``PressureRegulator`` links only: "controlling", "open" or "closed";
    # the other links at the controlled node; state changes so far.
    regulator_state: str = "controlling"
    regulator_capped: bool = False    # controlling, but limited by the wide-open valve this iteration
    regulator_cap_hits: int = 0       # iterations limited so in the current hydraulic pass
    regulator_can_control: bool = False   # wide open it overshot the set point: control is feasible
    regulator_links: tuple = ()
    regulator_switches: int = 0


@dataclass
class PumpState(PressureChangerState):
    component: Pump
    reynolds: float = 0.0


@dataclass
class HeatSourceState(PressureChangerState):
    component: HeatSource
    reynolds: float = 0.0
    temperature_c: float | None = None
    exit_temperature_c: float | None = None


@dataclass
class NetworkState:
    nodes: dict[int, NodeState]
    components: list[PressureChangerState]
    # separator node id → its ``SeparatorOutlet`` link states (set by the solvers)
    separators: dict = field(default_factory=dict)
    separator_splits: dict = field(default_factory=dict)   # node id → SeparatorSplit
    regulators: tuple = ()                                   # PressureRegulator link states
