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
    def _make_isothermal(network_state, T: float) -> tuple[str, ...]:
        """Hold every node, pipe and heater at ``T`` for an isothermal solve
        (``isothermal_temperature_c`` in the gas, black-oil and compositional
        settings): no energy equation, every node a fixed temperature, and
        ``Well`` tubing at the reservoir node's temperature (``T``).  Returns
        the warnings: heater/cooler duties and thermal boundaries at another
        temperature have no effect."""
        from angelica.core.components import Well, WellInflow
        from angelica.core.state import FittingState, PumpState

        warnings = []
        at = f"the solve is isothermal at {T:g} °C"
        for nid, node in network_state.nodes.items():
            if node.is_thermal_inlet and node.temperature_c is not None and abs(node.temperature_c - T) > 1e-9:
                warnings.append(f"node {nid}: thermal boundary {node.temperature_c:g} °C ignored, {at}")
            if node.thermal_gradient_dc_per_m:
                warnings.append(f"node {nid}: thermal gradient boundary ignored, {at}")
            node.temperature_c = T
            node.is_thermal_inlet = True
        pumps, valves, wells = [], [], []
        for link in network_state.components:
            cid = link.component.component_id
            if isinstance(link, HeatSourceState):
                if float(getattr(link.component, "power_w", 0.0) or 0.0) != 0.0:
                    warnings.append(f"{cid}: heater/cooler duty ignored, {at}")
                link.temperature_c = T
            elif isinstance(link, PipeState):
                link.temperature_c = T
            elif isinstance(link, PumpState):
                pumps.append(cid)
            elif isinstance(link.component, (Well, WellInflow)):
                wells.append(cid)
            elif isinstance(link, FittingState):
                valves.append(cid)
            link.isothermal = True     # read by the Well model (tubing at T)
        if pumps:
            warnings.append(f"{', '.join(pumps)}: heating by the pump/compressor not computed, {at}")
        if valves:
            warnings.append(f"{', '.join(valves)}: Joule-Thomson change across the valve(s) not computed, {at}")
        if wells:
            warnings.append(f"{', '.join(wells)}: tubing and inflow at {T:g} °C, no heat exchange with "
                            f"the formation (isothermal solve)")
        return tuple(warnings)

    def _energy_step(self, isothermal_c, network_state, fluid, scheme, T_ref: float,
                     pipe_enthalpy_terms: bool = True):
        """The energy solve, or the held temperatures of an isothermal solve."""
        if isothermal_c is not None:
            return self._held_temperatures(network_state, T_ref)
        from angelica.numerics.energy import solve_energy_system
        return solve_energy_system(network_state, fluid, scheme, T_ref=T_ref,
                                   pipe_enthalpy_terms=pipe_enthalpy_terms)

    @staticmethod
    def _held_temperatures(network_state, T_default: float):
        """Present node and pipe temperatures, in the form the energy solve
        returns them (pipe index counts pipes and heat sources)."""
        nodes = {nid: (n.temperature_c if n.temperature_c is not None else T_default)
                 for nid, n in network_state.nodes.items()}
        pipes = {}
        idx = 0
        for link in network_state.components:
            if isinstance(link, (PipeState, HeatSourceState)):
                if link.temperature_c is not None:
                    pipes[idx] = link.temperature_c
                idx += 1
        return nodes, pipes

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
            src = getattr(node, "source_mass_flow_kg_per_s", 0.0)
            mass_in += src                     # internal sources (lift gas) enter the network
            node_net_in[nid] += src
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
           real gas), plus the Joule-Thomson source of each pipe: the change
           of sensible ṁ·cp·T at constant enthalpy, negative when the gas cools.
        6. Ė_pot: −Σ ṁ·g·Δz along pipes (fluids with an enthalpy model).

        Balance: Ė_in + Q_src + W_shaft + Ė_JT + Ė_pot − Q_wall − Ė_out ≈ 0
        """

        from angelica.numerics.energy import (
            _has_enthalpy,
            _LinkAtTemperature,
            device_thermal_state,
            fluid_enthalpy,
        )

        real_basis = _has_enthalpy(fluid_model)

        def _h_rel(link, P, T: float) -> float:
            """ṁ-specific enthalpy above the reference (ideal gas at 0 °C)."""
            ref = (fluid_model.reference_enthalpy_for_link(link)
                   if hasattr(fluid_model, "reference_enthalpy_for_link") else 0.0)
            P_use = float(P) if P is not None else getattr(fluid_model, "reference_pressure_pa", 101_325.0)
            h = fluid_enthalpy(fluid_model, link, P_use, T)
            return (h if h is not None else 0.0) - ref

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
                if real_basis:
                    if leaves_node:
                        h_in_w += abs(flow) * _h_rel(link, node.pressure_pa, node.temperature_c)
                    else:
                        h_out_w += abs(flow) * _h_rel(link, node.pressure_pa,
                                                      _arriving_temperature(link, node.temperature_c))
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
        q_well_w = 0.0          # heat wells give to the formation
        e_pot_well_w = 0.0      # −ṁ·g·Δz of the well tubing
        for link_id, dt in devices.items():
            delta = abs(dt.mdot) * (dt.cp_out * dt.T_out - dt.cp_up * dt.T_up)
            if dt.kind == "well":
                q_well_w += dt.heat_loss_w
                e_pot_well_w += dt.potential_w
                if dt.lift_mdot > 0.0:          # lift gas: enters at the valve, leaves at the wellhead
                    h_in_w += dt.lift_mdot * dt.lift_h_supply_rel
                    if dt.down_node.is_boundary:
                        h_out_w += dt.lift_mdot * dt.lift_h_out_rel
            if dt.kind == "pump":
                if real_basis:
                    link = next(lk for lk in network_state.components if id(lk) == link_id)
                    delta = abs(dt.mdot) * (_h_rel(link, dt.down_node.pressure_pa, dt.T_out)
                                            - _h_rel(link, dt.up_node.pressure_pa, dt.T_up))
                w_shaft_w += delta
            elif dt.kind == "isenthalpic":
                e_jt_w += delta

        # Heat from link sources and pipe walls; pipe JT and potential-energy sources
        q_sources_w = 0.0
        q_wall_w = q_well_w
        e_pot_w = e_pot_well_w
        for link in network_state.components:
            if isinstance(link, HeatSourceState):
                q_sources_w += link.component.power_w
            elif isinstance(link, PipeState):
                e_jt_w += link.jt_source_w
                e_pot_w += link.gravity_source_w
                s_pipe = link.jt_source_w + link.gravity_source_w
                pipe = link.component
                U = pipe.u_overall_w_per_m2k
                if U <= 0.0:
                    continue
                # With a uniform source S the analytical profile relaxes to
                # T_amb + S/(U·π·D·L); the wall loss is U·π·D·L·LMTD about that
                # temperature plus S.
                ua = U * math.pi * pipe.diameter_m * pipe.length_m
                T_amb = pipe.ambient_temperature_c
                if pipe.n_thermal_segments == 1 and s_pipe != 0.0:
                    q_wall_w += s_pipe
                    T_amb = T_amb + s_pipe / ua
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
            potential_energy_kw=e_pot_w / 1000.0,
            enthalpy_basis="real" if real_basis else "sensible",
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
