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


def _rise_step(link, rise: float, f: float, P1: float, target: float, cap: float, limit_growth: bool) -> float:
    """Next pressure rise imposed by a compressor link.

    ``f`` = ln(P2/P1 at the nodes) − ln(target) increases with the imposed
    rise: a larger rise drives more flow, lowers P1 and raises P2.  The
    target is the specified ratio, (P1 + rise_spec)/P1, or the ratio the
    specified power delivers at the present flow.  One evaluation per outer
    iteration (each needs a hydraulic pass): a bracket on f is kept on the
    link and refined by false position with the Illinois modification.
    Without a bracket the step corrects the rise by the pressure it misses at
    the outlet (dP2/d rise ≈ 1); for a power specification the imposed ratio
    grows at most ×1.5 per step (zero flow asks for the ratio limit, whose
    flow emptied the suction line, v1.30.0 development).  The rise stays
    within [0, cap]."""
    lo, flo, hi, fhi, side = getattr(link, "rise_bracket", None) or (None, None, None, None, 0)
    if f < 0.0:
        lo, flo = rise, f
        if side == -1 and fhi is not None:
            fhi *= 0.5                     # Illinois: the other end retained twice
        side = -1
    else:
        hi, fhi = rise, f
        if side == 1 and flo is not None:
            flo *= 0.5
        side = 1
    if lo is not None and hi is not None and (hi - lo) <= 1e-6 * max(abs(hi), 1e3):
        # Collapsed onto one rise with f of both signs, or inverted (f < 0 at a
        # larger rise than f > 0): f moved between passes (temperatures, or
        # the flow through other links, shift the pressures); the old end is
        # stale.  (Inverted, the step was clamped to the same rise for ever —
        # a flow-specified inlet, found in review.)
        lo = flo = hi = fhi = None
        side = 0
    link.rise_bracket = (lo, flo, hi, fhi, side)
    if lo is not None and hi is not None and flo is not None and fhi is not None:
        x = lo - flo * (hi - lo) / (fhi - flo)
        return min(max(x, lo), hi, cap)
    new = rise - math.expm1(f) * P1 * target           # outlet pressure missed: P2 − target·P1
    if limit_growth:
        new = min(new, 1.5 * (P1 + rise) - P1)
    return min(max(new, 0.0), cap)


# Lowest node pressure of the solvers that work in absolute pressure (gas,
# black-oil, compositional): the equation of state divides by P.
ABSOLUTE_PRESSURE_FLOOR_PA = 1.0e3

# Relative change of a power-specified compressor's pressure rise below which
# the outer loop counts it as settled.
COMPRESSOR_RISE_TOLERANCE = 1e-4


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
                if float(getattr(link.component, "power_w", 0.0) or 0.0) != 0.0 \
                        or link.component.outlet_temperature_c is not None:
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
    def _set_compressor_rises(network_state, fluid, T_default: float = 20.0) -> float:
        """Pressure rise each ``Compressor`` imposes, adjusted until the node
        pressures meet its specification (at the present flow and inlet
        temperature for a power); returns the largest mismatch
        |ln(P2/P1) − ln(target)| at the nodes."""
        from angelica.closures.compressor import CompressorFluid, outlet_pressure_for_power
        from angelica.core.state import CompressorState
        change = 0.0
        for link in network_state.components:
            if not isinstance(link, CompressorState):
                continue
            comp = link.component
            P1 = float(link.start_node.pressure_pa or 0.0)
            P2 = float(link.end_node.pressure_pa or 0.0)
            if P1 <= 0.0 or P2 <= 0.0:
                continue
            if comp.specification == "pressure ratio":
                target = comp.pressure_ratio
            elif comp.specification == "pressure rise":
                target = (P1 + comp.pressure_rise_pa) / P1
            else:
                T1 = link.start_node.temperature_c if link.start_node.temperature_c is not None else T_default
                target = outlet_pressure_for_power(CompressorFluid(fluid, link), comp, P1, T1,
                                                   max(float(link.mass_flow_kg_per_s), 0.0)) / P1
            target = max(target, 1.0)
            if getattr(link, "rise_bracket", None) is None and link.compressor_rise_pa == 0.0 \
                    and comp.specification != "power":
                link.compressor_rise_pa = (target - 1.0) * P1      # first guess: a transparent link
            f = math.log(P2 / P1) - math.log(target)
            # The outer loop is converged only when the pressures at the nodes
            # meet the specification, not when the imposed rise stops moving —
            # except for a machine the network pushes beyond its specification
            # (f > 0) with no rise left and no flow: it is idle, a settled
            # state (e.g. a ratio-1.8 unit beside a ratio-2 one in parallel, or
            # an outlet above what it can reach).
            idle = f > 0.0 and link.compressor_rise_pa <= 0.0 and float(link.mass_flow_kg_per_s) <= 1e-9
            if not idle:
                change = max(change, abs(f))
            cap = (max(comp.max_pressure_ratio, target) * 1.5 - 1.0) * P1
            link.compressor_rise_pa = _rise_step(link, max(link.compressor_rise_pa, 0.0), f, P1, target, cap,
                                                 comp.specification == "power")
        return change

    @staticmethod
    def _compressor_warnings(network_state, fluid, T_default: float = 20.0) -> tuple[str, ...]:
        """Compressors fed with a stream that is not (almost) all gas."""
        from angelica.closures.compressor import CompressorFluid
        from angelica.core.state import CompressorState
        out = []
        for link in network_state.components:
            if not isinstance(link, CompressorState) or float(link.mass_flow_kg_per_s) <= 1e-12:
                continue
            P1 = float(link.start_node.pressure_pa or 0.0)
            T1 = link.start_node.temperature_c if link.start_node.temperature_c is not None else T_default
            try:
                note = CompressorFluid(fluid, link).liquid_at(P1, T1)
            except Exception:       # a check must never break a converged solve
                note = None
            if note:
                out.append(f"{link.component.component_id}: the feed is not all gas at suction ({note}); "
                           "the compressor model is a single-phase gas's")
        return tuple(out)

    @staticmethod
    def _compressor_results(network_state, fluid, T_default: float = 20.0) -> tuple:
        """Power, pressure ratio and temperatures of every ``Compressor``."""
        from angelica.closures.compressor import CompressorFluid, compress
        from angelica.core.results import CompressorResult
        from angelica.core.state import CompressorState
        out = []
        for link in network_state.components:
            if not isinstance(link, CompressorState):
                continue
            comp = link.component
            m = float(link.mass_flow_kg_per_s)
            P1, P2 = float(link.start_node.pressure_pa or 0.0), float(link.end_node.pressure_pa or 0.0)
            T1 = link.start_node.temperature_c if link.start_node.temperature_c is not None else T_default
            if m <= 1e-12 or P2 <= P1 or P1 <= 0.0:
                state = "no flow" if m <= 1e-12 else "no pressure rise"
                out.append(CompressorResult(comp.component_id, comp.specification, state, m, P1, P2,
                                            P2 / P1 if P1 > 0 else float("nan"), T1, T1, T1, 0.0, 0.0))
                continue
            res = compress(CompressorFluid(fluid, link), P1, T1, P2, comp.efficiency, comp.efficiency_type,
                           comp.n_steps)
            state = "running"
            if comp.specification == "power" and P2 / P1 >= comp.max_pressure_ratio * (1.0 - 1e-6):
                state = "at max pressure ratio"
            out.append(CompressorResult(comp.component_id, comp.specification, state, m, P1, P2, P2 / P1, T1,
                                        res.t_out_c, res.t_out_isentropic_c, m * res.actual_work_j_per_kg,
                                        m * res.isentropic_work_j_per_kg))
        return tuple(out)

    @staticmethod
    def _heater_results(network_state, isothermal: bool = False) -> tuple:
        """Duty and outlet temperature of every ``HeatSource`` (heater/cooler)."""
        from angelica.core.results import HeaterResult
        out = []
        for link in network_state.components:
            if not isinstance(link, HeatSourceState):
                continue
            comp = link.component
            m = float(link.mass_flow_kg_per_s)
            up, down = (link.start_node, link.end_node) if m >= 0.0 else (link.end_node, link.start_node)
            T_in = up.temperature_c if up.temperature_c is not None else float("nan")
            T_out = (link.exit_temperature_c if link.exit_temperature_c is not None and not isothermal
                     else T_in)
            if isothermal or abs(m) < 1e-30 or link.duty_w is None:
                duty, state = 0.0, "no flow" if abs(m) < 1e-30 else "isothermal (duty ignored)"
            else:
                duty = float(link.duty_w)
                state = ("fixed duty" if comp.outlet_temperature_c is None
                         else "duty limited" if link.duty_limited else "at outlet temperature")
            out.append(HeaterResult(component_id=comp.component_id, state=state, duty_w=duty,
                                    inlet_temperature_c=float(T_in), outlet_temperature_c=float(T_out),
                                    outlet_temperature_spec_c=comp.outlet_temperature_c, mass_flow_kg_per_s=m))
        return tuple(out)

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
        2. Q_src: Σ duty of HeatSource links (fixed, or from the outlet
           temperature set point — as set by the energy solve).
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
                q_sources_w += link.duty_w if link.duty_w is not None else link.component.power_w
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
