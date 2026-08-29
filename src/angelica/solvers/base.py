from __future__ import annotations

from abc import ABC, abstractmethod
from typing import TYPE_CHECKING, Callable, Optional

if TYPE_CHECKING:
    from angelica.core.case import NetworkCase
    from angelica.core.results import SolveResult

import math

from angelica.core.results import GlobalBalance, GlobalEnergyBalance
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
    def _compute_global_balance(network_state) -> GlobalBalance:
        # Net flow INTO each node (positive = fluid arriving at this node)
        node_net_in: dict[int, float] = {nid: 0.0 for nid in network_state.nodes}
        for link in network_state.components:
            flow = link.mass_flow_kg_per_s
            node_net_in[link.start_node.node_id] -= flow
            node_net_in[link.end_node.node_id] += flow

        mass_in = 0.0
        mass_out = 0.0
        for nid, node in network_state.nodes.items():
            if not node.is_boundary:
                continue
            net = node_net_in[nid]
            if net < 0:
                mass_in += -net   # fluid injected into the network
            else:
                mass_out += net   # fluid absorbed from the network

        return GlobalBalance(
            mass_inlet_kg_per_s=mass_in,
            mass_outlet_kg_per_s=mass_out,
        )

    @staticmethod
    def _compute_global_energy_balance(
        network_state,
        fluid_model,
        T_ref: float = 20.0,
    ) -> GlobalEnergyBalance:
        """Compute the steady-state global energy balance (first law).

        1. Enthalpy entering via inlet boundary nodes: Ė_in = Σ ṁ·cp·T
        2. Enthalpy leaving via outlet boundary nodes: Ė_out = Σ ṁ·cp·T
        3. Heat added by HeatSource links: Q_src = Σ power_w
        4. Heat lost through pipe walls: Q_wall = Σ U·π·D·L·(T_pipe − T_amb)

        Balance: Ė_in + Q_src − Q_wall − Ė_out ≈ 0
        """

        class _TempProxy:
            __slots__ = ("temperature_c",)
            def __init__(self, t: float) -> None:
                self.temperature_c = t

        # Net flow into each node
        node_net_in: dict[int, float] = {nid: 0.0 for nid in network_state.nodes}
        for link in network_state.components:
            flow = link.mass_flow_kg_per_s
            node_net_in[link.start_node.node_id] -= flow
            node_net_in[link.end_node.node_id] += flow

        # Boundary enthalpy fluxes
        h_in_w = 0.0
        h_out_w = 0.0
        for nid, node in network_state.nodes.items():
            if not node.is_boundary:
                continue
            T_node = node.temperature_c
            if T_node is None:
                continue
            net = node_net_in[nid]
            if abs(net) < 1e-30:
                continue
            try:
                cp = fluid_model.specific_heat_for_link(_TempProxy(T_node))
            except NotImplementedError:
                continue
            rate = abs(net) * cp * T_node  # W  (0 °C reference)
            if net < 0:           # fluid injected into network
                h_in_w += rate
            else:                 # fluid absorbed from network
                h_out_w += rate

        # Heat from link sources and pipe walls
        q_sources_w = 0.0
        q_wall_w = 0.0
        for link in network_state.components:
            if isinstance(link, HeatSourceState):
                q_sources_w += link.component.power_w
            elif isinstance(link, PipeState):
                pipe = link.component
                U = pipe.heat_transfer_coefficient_w_per_m2k
                if U <= 0.0:
                    continue
                T_amb = pipe.ambient_temperature_c
                if pipe.n_thermal_segments == 1:
                    # NTU analytical bypass: link.temperature_c IS the exit temperature.
                    # Use LMTD: Q = U·π·D·L·LMTD, exact for the NTU bypass path.
                    mdot = abs(link.mass_flow_kg_per_s)
                    if mdot < 1e-30:
                        continue
                    if link.mass_flow_kg_per_s >= 0:
                        _t = network_state.nodes[link.start_node.node_id].temperature_c
                        T_up = _t if _t is not None else T_ref
                    else:
                        _t = network_state.nodes[link.end_node.node_id].temperature_c
                        T_up = _t if _t is not None else T_ref
                    T_exit = link.temperature_c if link.temperature_c is not None else T_ref
                    dT1 = T_up - T_amb
                    dT2 = T_exit - T_amb
                    if abs(dT1) < 1e-12 or abs(dT2) < 1e-12 or dT1 * dT2 <= 0.0:
                        # Near-ambient or sign flip: fall back to enthalpy change
                        try:
                            cp = fluid_model.specific_heat_for_link(_TempProxy(T_up))
                        except NotImplementedError:
                            continue
                        q_wall_w += mdot * cp * (T_up - T_exit)
                    else:
                        ratio = dT1 / dT2
                        lmtd = (dT1 - dT2) / math.log(ratio) if abs(ratio - 1.0) > 1e-9 else dT1
                        q_wall_w += U * math.pi * pipe.diameter_m * pipe.length_m * lmtd
                else:
                    # Multi-segment FV: source terms act on n_internal = n_segs−1 nodes,
                    # each covering dx = L/n_segs. Effective wall length = L·n_internal/n_segs.
                    n_segs = max(pipe.n_thermal_segments, 2)
                    n_internal = n_segs - 1
                    fv_length = pipe.length_m * n_internal / n_segs
                    T_pipe = link.temperature_c if link.temperature_c is not None else T_ref
                    q_wall_w += U * math.pi * pipe.diameter_m * fv_length * (T_pipe - T_amb)

        return GlobalEnergyBalance(
            enthalpy_in_kw=h_in_w / 1000.0,
            enthalpy_out_kw=h_out_w / 1000.0,
            heat_sources_kw=q_sources_w / 1000.0,
            heat_wall_loss_kw=q_wall_w / 1000.0,
        )

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
                    _ts = network_state.nodes[component_state.start_node.node_id].temperature_c
                    _te = network_state.nodes[component_state.end_node.node_id].temperature_c
                    start_temp = _ts if _ts is not None else default_temperature_c
                    end_temp   = _te if _te is not None else default_temperature_c
                    component_state.temperature_c = 0.5 * (start_temp + end_temp)
                pipe_idx += 1
