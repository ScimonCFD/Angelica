from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Callable, Optional

if TYPE_CHECKING:
    from angelica.core.case import NetworkCase
    from angelica.core.results import SolveResult

from angelica.core.state import HeatSourceState, PipeState


class BaseSolver(ABC):
    @abstractmethod
    def solve(
        self,
        case: NetworkCase,
        progress_callback: Optional[Callable] = None,
    ) -> SolveResult:
        raise NotImplementedError

    @staticmethod
    def _require_fixed_temperature(case: NetworkCase) -> None:
        """Raise if no fixed-temperature thermal boundary is defined.

        The energy system is singular without at least one Dirichlet node.
        """
        for tb in case.thermal_inlets:
            if tb.bc_type == "fixed_temperature":
                return
        raise ValueError(
            "The thermal solver requires at least one ThermalBoundary with "
            "bc_type='fixed_temperature'. Without it the energy system has no "
            "unique solution (pure Neumann problem)."
        )

    @staticmethod
    def _initial_temperature(case: NetworkCase) -> float:
        for tb in case.thermal_inlets:
            if tb.bc_type == "fixed_temperature":
                return tb.temperature_c
        for comp in case.components:
            if hasattr(comp, "ambient_temperature_c"):
                return comp.ambient_temperature_c
        return 20.0

    @staticmethod
    def _update_component_temperatures(
        network_state,
        default_temperature_c: float,
        pipe_mean_temps: dict[int, float] | None = None,
    ) -> None:
        pipe_idx = 0
        for component_state in network_state.components:
            if isinstance(component_state, (PipeState, HeatSourceState)):
                if pipe_mean_temps is not None and pipe_idx in pipe_mean_temps:
                    component_state.temperature_c = pipe_mean_temps[pipe_idx]
                else:
                    start_temp = (
                        network_state.nodes[component_state.start_node.node_id].temperature_c
                        or default_temperature_c
                    )
                    end_temp = (
                        network_state.nodes[component_state.end_node.node_id].temperature_c
                        or default_temperature_c
                    )
                    component_state.temperature_c = 0.5 * (start_temp + end_temp)
                pipe_idx += 1
