from __future__ import annotations

from dataclasses import dataclass, field

from .components import PressureChanger
from angelica.properties.base import FluidModel


@dataclass(frozen=True)
class InletFluidBC:
    """Black-oil fluid composition at a pressure or flow inlet node.

    Specifies the four characterisation parameters that define the produced
    fluid from a particular reservoir or well.  When multiple inlets with
    different compositions feed a network, the solver propagates each inlet's
    fluid through the network and mixes at junctions by mass-weighted average.

    Args:
        node_id: ID of the inlet node this composition is assigned to.
        api_gravity: Stock-tank oil API gravity (°API).
        gas_gravity: Gas specific gravity relative to air (–).
        gor_sc_m3_per_m3: Gas-oil ratio at standard conditions (m³/m³).
        wor_sc_m3_per_m3: Water-oil ratio at standard conditions (m³/m³).
    """
    node_id: int
    api_gravity: float
    gas_gravity: float
    gor_sc_m3_per_m3: float
    wor_sc_m3_per_m3: float


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
    """Thermal boundary condition for a node in a non-isothermal network.

    Args:
        node_id: ID of the node where this boundary condition is applied.
        temperature_c: Prescribed temperature (°C). Used only when
            ``bc_type="fixed_temperature"``.
        bc_type: Type of thermal boundary condition:

            * ``"fixed_temperature"`` — Dirichlet BC. Sets the nodal
              temperature to ``temperature_c`` and holds it fixed throughout
              the outer iteration.  Use this at supply inlets where the fluid
              enters at a known temperature.

            * ``"zero_gradient"`` — Neumann BC with zero normal gradient
              (∂T/∂n = 0). The solver imputes the outlet temperature from the
              upwind FV solution without imposing a value.  Use this at
              pressure outlets when the exit temperature is unknown.

            * ``"fixed_gradient"`` — Neumann BC with a prescribed non-zero
              gradient ``gradient_dc_per_m`` [°C/m].  Rarely needed; defaults
              to zero which is equivalent to ``"zero_gradient"``.

        gradient_dc_per_m: Temperature gradient [°C/m] for
            ``"fixed_gradient"`` BCs.  Ignored for other bc types.
    """

    _VALID_BC_TYPES = ("fixed_temperature", "zero_gradient", "fixed_gradient")

    node_id: int
    temperature_c: float = 0.0
    bc_type: str = "fixed_temperature"
    gradient_dc_per_m: float = 0.0

    def __post_init__(self) -> None:
        if self.bc_type not in self._VALID_BC_TYPES:
            raise ValueError(
                f"ThermalBoundary bc_type must be one of {self._VALID_BC_TYPES!r}, "
                f"got {self.bc_type!r}"
            )


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
    inlet_fluid_bcs: tuple[InletFluidBC, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        _tuple_fields = (
            "pressure_inlets", "pressure_outlets", "components",
            "flow_inlets", "flow_outlets", "node_ids", "thermal_inlets",
            "inlet_fluid_bcs",
        )
        for attr in _tuple_fields:
            val = getattr(self, attr)
            if not isinstance(val, tuple):
                object.__setattr__(self, attr, tuple(val))

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
