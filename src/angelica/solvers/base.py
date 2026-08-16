from __future__ import annotations

from abc import ABC, abstractmethod

from angelica.core.state import HeatSourceState, PipeState


class BaseSolver(ABC):
    @abstractmethod
    def solve(self, case):
        raise NotImplementedError

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
