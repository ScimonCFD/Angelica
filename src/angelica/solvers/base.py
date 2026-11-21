from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Callable
from typing import TYPE_CHECKING

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
        progress_callback: Callable | None = None,
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
        node_net_in: dict[int, float] = dict.fromkeys(network_state.nodes, 0.0)
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

        Quantities are the ones the energy equation conserves, ṁ·cp·T with the
        heat capacity of each link (0 °C reference):

        1. Ė_in / Ė_out: Σ ṁ·cp_link·T_node over the links attached to inlet /
           outlet boundary nodes (cp of a pipe at its mean temperature, of a
           device at its inlet / outlet state — as in the assembly).
        2. Q_src: Σ power_w of HeatSource links.
        3. Q_wall: Σ U·π·D·L·(T_pipe − T_amb), computed independently of the
           assembly so that discretisation errors in the wall term show up.
        4. W_shaft: Σ ṁ·(cp·T)_out − ṁ·(cp·T)_in over pumps/compressors
           (= ṁ·ΔP/ρ for liquids).
        5. Ė_JT: the same difference over isenthalpic devices (valves in a
           real gas): the change of sensible ṁ·cp·T at constant enthalpy,
           negative when the gas cools.

        Balance: Ė_in + Q_src + W_shaft + Ė_JT − Q_wall − Ė_out ≈ 0
        """

        from angelica.numerics.energy import _LinkAtTemperature, device_thermal_state

        devices = {}
        for link in network_state.components:
            dt = device_thermal_state(link, fluid_model, T_ref)
            if dt is not None:
                devices[id(link)] = dt

        def _cp_at_node(link, node_id: int, T_node: float) -> float:
            dt = devices.get(id(link))
            if dt is not None:
                return dt.cp_up if dt.up_node.node_id == node_id else dt.cp_out
            T_link = link.temperature_c if link.temperature_c is not None else T_node
            return fluid_model.specific_heat_for_link(_LinkAtTemperature(link, T_link))

        def _arriving_temperature(link, T_node: float) -> float:
            """Temperature of the fluid a link delivers into a boundary node.

            Usually equal to the node temperature, but not when the node has a
            fixed temperature and receives back-flow (e.g. a well taking flow
            from the manifold): the fluid arrives at the pipe's exit temperature.
            """
            dt = devices.get(id(link))
            if dt is not None:
                return dt.T_out
            t_exit = getattr(link, "exit_temperature_c", None)
            return t_exit if t_exit is not None else T_node

        # Boundary enthalpy fluxes, link by link
        h_in_w = 0.0
        h_out_w = 0.0
        for link in network_state.components:
            flow = float(link.mass_flow_kg_per_s)
            if abs(flow) < 1e-30:
                continue
            for node, leaves_node in ((link.start_node, flow > 0.0), (link.end_node, flow < 0.0)):
                if not node.is_boundary or node.temperature_c is None:
                    continue
                try:
                    cp = _cp_at_node(link, node.node_id, node.temperature_c)
                except NotImplementedError:
                    continue
                if leaves_node:      # fluid injected into the network here
                    h_in_w += abs(flow) * cp * node.temperature_c      # W (0 °C reference)
                else:                # fluid absorbed by the boundary here
                    h_out_w += abs(flow) * cp * _arriving_temperature(link, node.temperature_c)

        # Shaft work and Joule-Thomson term of devices
        w_shaft_w = 0.0
        e_jt_w = 0.0
        for dt in devices.values():
            delta = abs(dt.mdot) * (dt.cp_out * dt.T_out - dt.cp_up * dt.T_up)
            if dt.kind == "pump":
                w_shaft_w += delta
            elif dt.kind == "isenthalpic":
                e_jt_w += delta

        # Heat from link sources and pipe walls
        q_sources_w = 0.0
        q_wall_w = 0.0
        for link in network_state.components:
            if isinstance(link, HeatSourceState):
                q_sources_w += link.component.power_w
            elif isinstance(link, PipeState):
                pipe = link.component
                U = pipe.u_overall_w_per_m2k
                if U <= 0.0:
                    continue
                T_amb = pipe.ambient_temperature_c
                if pipe.n_thermal_segments == 1:
                    # NTU analytical bypass: link.exit_temperature_c is the pipe exit temperature.
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
                    T_exit = (link.exit_temperature_c if link.exit_temperature_c is not None
                              else link.temperature_c if link.temperature_c is not None else T_ref)
                    dT1 = T_up - T_amb
                    dT2 = T_exit - T_amb
                    if abs(dT1) < 1e-12 or abs(dT2) < 1e-12 or dT1 * dT2 <= 0.0:
                        # Near-ambient or sign flip: fall back to enthalpy change
                        try:
                            cp = fluid_model.specific_heat_for_link(_LinkAtTemperature(link, T_up))
                        except NotImplementedError:
                            continue
                        q_wall_w += mdot * cp * (T_up - T_exit)
                    else:
                        ratio = dT1 / dT2
                        lmtd = (dT1 - dT2) / math.log(ratio) if abs(ratio - 1.0) > 1e-9 else dT1
                        q_wall_w += U * math.pi * pipe.diameter_m * pipe.length_m * lmtd
                else:
                    # Multi-segment FV: the full wall area U·π·D·L is shared by the
                    # n_internal nodes; link.temperature_c is their mean temperature.
                    T_pipe = link.temperature_c if link.temperature_c is not None else T_ref
                    q_wall_w += U * math.pi * pipe.diameter_m * pipe.length_m * (T_pipe - T_amb)

        return GlobalEnergyBalance(
            enthalpy_in_kw=h_in_w / 1000.0,
            enthalpy_out_kw=h_out_w / 1000.0,
            heat_sources_kw=q_sources_w / 1000.0,
            heat_wall_loss_kw=q_wall_w / 1000.0,
            shaft_work_kw=w_shaft_w / 1000.0,
            joule_thomson_kw=e_jt_w / 1000.0,
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
