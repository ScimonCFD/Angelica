from __future__ import annotations

import math
import warnings
from dataclasses import dataclass
from typing import Any

import numpy as np

from angelica.closures.convection_scheme import ConvectionScheme
from angelica.closures.gravity import GRAVITY_M_PER_S2
from angelica.core.state import FittingState, HeatSourceState, NetworkState, PipeState, PumpState


def _load_scipy():
    try:
        # Some local Python environments emit a conservative compatibility
        # warning during SciPy import even though the sparse solve path used by
        # Angelica works correctly. Keep tutorial and GUI output clean, while
        # still failing hard on a genuine import error.
        with warnings.catch_warnings():
            warnings.filterwarnings(
                "ignore",
                message=r"A NumPy version .* is required for this version of SciPy",
                category=UserWarning,
                module="scipy",
            )
            import scipy.sparse as sp
            import scipy.sparse.linalg as spla
    except ImportError as exc:
        raise RuntimeError(
            "The non-isothermal solver requires SciPy, which is not installed. "
            "Install it with: pip install scipy"
        ) from exc
    return sp, spla


def fluid_enthalpy(fluid_model, link_state, pressure_pa: float, temperature_c: float) -> float | None:
    """Specific enthalpy (J/kg) of the fluid in a link, or None when the fluid
    model has no enthalpy (incompressible liquids, black-oil)."""
    if hasattr(fluid_model, "enthalpy_for_link"):
        return fluid_model.enthalpy_for_link(link_state, pressure_pa, temperature_c)
    if hasattr(fluid_model, "enthalpy_j_per_kg"):
        return fluid_model.enthalpy_j_per_kg(pressure_pa, temperature_c)
    return None


def pipe_enthalpy_sources(ps, fluid_model, cp: float, T_ref: float = 20.0) -> tuple[float, float]:
    """(Joule-Thomson, potential-energy) sources of a pipe in W, + heats the fluid.

    The temperature equation transports ṁ·cp·T with the link heat capacity
    ``cp``; the first law along the pipe is ṁ·Δh = −Q_wall − ṁ·g·Δz.  The part
    of Δh that ṁ·cp·ΔT does not carry is returned as a source:

        S_JT  = −ṁ·{[h(P_out, T_out) − h(P_in, T_in)] − cp·(T_out − T_in)}
        S_pot = −ṁ·g·Δz                                  (flow direction)

    For a real gas S_JT is the Joule-Thomson cooling (friction heating is
    included through the actual outlet pressure) together with the real-gas
    heat-capacity departure; for a mixture with phase change it also carries
    latent heat.  T_out is the pipe's exit temperature from the previous
    energy solve (the inlet temperature on the first), so at convergence the
    first law holds exactly.  Zero for fluids without an enthalpy model.
    """
    if not isinstance(ps, PipeState):
        return 0.0, 0.0
    mdot = abs(float(ps.mass_flow_kg_per_s))
    if mdot < 1e-30:
        return 0.0, 0.0
    up, down = (ps.start_node, ps.end_node) if ps.mass_flow_kg_per_s >= 0.0 else (ps.end_node, ps.start_node)
    P_up, P_down = up.pressure_pa, down.pressure_pa
    if P_up is None or P_down is None:
        return 0.0, 0.0
    T_in = up.temperature_c if up.temperature_c is not None else T_ref
    T_out = ps.exit_temperature_c if ps.exit_temperature_c is not None else T_in
    h_in = fluid_enthalpy(fluid_model, ps, float(P_up), T_in)
    if h_in is None:
        return 0.0, 0.0
    h_out = fluid_enthalpy(fluid_model, ps, float(P_down), T_out)
    if h_out is None:
        return 0.0, 0.0
    dz = ps.component.height_change_m * (1.0 if ps.mass_flow_kg_per_s >= 0.0 else -1.0)
    return -mdot * ((h_out - h_in) - cp * (T_out - T_in)), -mdot * GRAVITY_M_PER_S2 * dz


def solve_energy_system(
    network_state: NetworkState,
    fluid_model,
    convection_scheme: ConvectionScheme,
    T_ref: float = 20.0,
    pipe_enthalpy_terms: bool = True,
    max_source_iterations: int = 30,
    source_tolerance_k: float = 2e-3,
) -> tuple[dict[int, float], dict[int, float]]:
    """Steady energy equation, with the pipe enthalpy sources made consistent.

    The Joule-Thomson / latent-heat source of a pipe depends on its exit
    temperature.  With the flow field fixed, the energy system is re-solved
    with updated sources until node and exit temperatures change by less than
    ``source_tolerance_k`` — cheap compared with a hydraulic solve, and it
    keeps the outer (hydraulics ↔ temperature) loop from crawling.  Node
    temperatures on the network state are restored afterwards; the caller
    applies (and relaxes) the returned ones.
    """
    saved = {nid: node.temperature_c for nid, node in network_state.nodes.items()}
    try:
        result = _solve_energy_once(network_state, fluid_model, convection_scheme, T_ref, pipe_enthalpy_terms)
        has_sources = pipe_enthalpy_terms and any(
            isinstance(ps, PipeState) and (ps.jt_source_w != 0.0 or ps.gravity_source_w != 0.0)
            for ps in network_state.components
        )
        if not has_sources or not _has_enthalpy(fluid_model):
            return result
        for _ in range(max_source_iterations - 1):
            node_temps, _means = result
            exits_before = [getattr(ps, "exit_temperature_c", None) for ps in network_state.components]
            change = 0.0
            for nid, T_new in node_temps.items():
                node = network_state.nodes[nid]
                if not node.is_thermal_inlet:
                    if node.temperature_c is not None:
                        change = max(change, abs(T_new - node.temperature_c))
                    node.temperature_c = T_new
            result = _solve_energy_once(network_state, fluid_model, convection_scheme, T_ref, pipe_enthalpy_terms)
            for ps, before in zip(network_state.components, exits_before):
                after = getattr(ps, "exit_temperature_c", None)
                if before is not None and after is not None:
                    change = max(change, abs(after - before))
            if change < source_tolerance_k:
                break
        return result
    finally:
        for nid, T in saved.items():
            network_state.nodes[nid].temperature_c = T


_SOURCE_RELAXATION = 0.5


def _has_enthalpy(fluid_model) -> bool:
    return hasattr(fluid_model, "enthalpy_for_link") or hasattr(fluid_model, "enthalpy_j_per_kg")


def _solve_energy_once(
    network_state: NetworkState,
    fluid_model,
    convection_scheme: ConvectionScheme,
    T_ref: float = 20.0,
    pipe_enthalpy_terms: bool = True,
) -> tuple[dict[int, float], dict[int, float]]:
    """Solve the steady-state energy equation for the pipe network.

    Returns:
        (node_temperatures, component_mean_temperatures) where:
        - node_temperatures: maps node_id → temperature_c for all network nodes.
        - component_mean_temperatures: maps pipe_states list index → mean temperature_c
          computed from the internal FV nodes, used for fluid property evaluation.

    Thermal BC types (from NodeState):
      - is_thermal_inlet = True       → Dirichlet: T = T_prescribed
      - thermal_gradient_dc_per_m = g → Neumann: dT/dx = g at boundary face
      - neither                        → junction mixing equation (≈ zero gradient)

    For fluids with an enthalpy model (compressible gas, compositional EOS)
    and ``pipe_enthalpy_terms=True``, pipes carry the Joule-Thomson and
    potential-energy sources of ``pipe_enthalpy_sources`` (stored on the pipe
    state for the global energy balance).

    The Neumann BC modifies the pipe FV at the boundary face:
      - removes diffusive coupling to the junction node
      - adds the prescribed gradient flux to the RHS
    The junction equation for a Neumann node extrapolates from the first
    interior node: T_junction = T_interior ± g * dx.
    """
    sp, spla = _load_scipy()

    sorted_node_ids = sorted(network_state.nodes.keys())
    N_nodes = len(sorted_node_ids)
    node_index = {nid: i for i, nid in enumerate(sorted_node_ids)}

    # ── index internal pipe/heat-source nodes ─────────────────────────────
    pipe_states = [ls for ls in network_state.components if isinstance(ls, (PipeState, HeatSourceState))]

    pipe_internal_offset: list[int] = []
    total_internal = 0
    for ps in pipe_states:
        pipe_internal_offset.append(N_nodes + total_internal)
        n_segs = max(ps.component.n_thermal_segments, 2)
        total_internal += n_segs - 1

    N_total = N_nodes + total_internal

    # ── collect Neumann gradient nodes ────────────────────────────────────
    # gradient_node_map: nid -> prescribed dT/dx (°C/m) at boundary face
    gradient_node_map: dict[int, float] = {}
    for node in network_state.nodes.values():
        if node.thermal_gradient_dc_per_m is not None and not node.is_thermal_inlet:
            gradient_node_map[node.node_id] = node.thermal_gradient_dc_per_m

    # For each gradient node, find the first connected pipe with n_internal > 0
    # so we can write a junction extrapolation equation.
    # Value: (pipe_idx, is_start_node, dx, offset, n_internal)
    gradient_node_pipe: dict[int, tuple[int, bool, float, int, int]] = {}
    for pipe_idx, ps in enumerate(pipe_states):
        n_segs_p = max(ps.component.n_thermal_segments, 2)
        n_internal_p = n_segs_p - 1
        if n_internal_p == 0:
            continue
        dx_p = ps.component.length_m / n_segs_p
        offset_p = pipe_internal_offset[pipe_idx]
        for nid, is_start in (
            (ps.start_node.node_id, True),
            (ps.end_node.node_id, False),
        ):
            if nid in gradient_node_map and nid not in gradient_node_pipe:
                gradient_node_pipe[nid] = (pipe_idx, is_start, dx_p, offset_p, n_internal_p)

    # ── build sparse matrix and RHS ───────────────────────────────────────
    rows: list[int] = []
    cols: list[int] = []
    vals: list[float] = []
    rhs = np.zeros(N_total)

    def add(r: int, c: int, v: float) -> None:
        rows.append(r)
        cols.append(c)
        vals.append(v)

    # NTU of pipes solved with the analytical single-segment path, used to
    # report their length-averaged temperature.
    ntu_paths: dict[int, tuple[float, float, float]] = {}   # pipe_idx → (NTU, T_amb_eff, S/|ṁcp|)

    # ── FV equations for pipe internal nodes ──────────────────────────────
    for pipe_idx, ps in enumerate(pipe_states):
        n_segs = max(ps.component.n_thermal_segments, 2)
        n_internal = n_segs - 1
        offset = pipe_internal_offset[pipe_idx]

        pipe = ps.component
        L = pipe.length_m
        dx = L / n_segs
        A = ps.area_m2
        D_pipe = pipe.diameter_m

        if isinstance(ps, HeatSourceState):
            U = 0.0
            T_amb = T_ref  # U=0 so this is unused; kept for structural symmetry with pipe branch
        else:
            assert isinstance(ps, PipeState)
            U = ps.component.u_overall_w_per_m2k
            T_amb = ps.component.ambient_temperature_c

        T_repr = ps.temperature_c if ps.temperature_c is not None else T_ref
        _ls = _LinkAtTemperature(ps, T_repr)

        cp = fluid_model.specific_heat_for_link(_ls)
        k_fl = fluid_model.thermal_conductivity_for_link(_ls)
        mdot = float(ps.mass_flow_kg_per_s)

        F = mdot * cp          # convective strength (W/K), signed
        D = k_fl * A / dx      # diffusive conductance (W/K)

        # Joule-Thomson and potential-energy sources (W), spread along the pipe
        if pipe_enthalpy_terms and isinstance(ps, PipeState):
            jt_new, ps.gravity_source_w = pipe_enthalpy_sources(ps, fluid_model, cp, T_ref)
            # Under-relaxed: with phase change the source falls steeply as the
            # exit temperature rises (latent heat), and the plain fixed point
            # alternates around the solution.
            ps.jt_source_w += _SOURCE_RELAXATION * (jt_new - ps.jt_source_w)
        elif isinstance(ps, PipeState):
            ps.jt_source_w = ps.gravity_source_w = 0.0
        S_pipe = ps.jt_source_w + ps.gravity_source_w if isinstance(ps, PipeState) else 0.0

        # ── Analytical NTU path for single-segment pipes ───────────────────
        # For pipes with n_thermal_segments=1 the FV discretisation places the
        # sole interior node at the midpoint, causing exit temperatures to
        # converge to exp(-NTU/2) instead of exp(-NTU). We substitute the
        # exact NTU solution: T_exit = decay·T_in + (1-decay)·T_amb.
        if (not isinstance(ps, HeatSourceState)
                and ps.component.n_thermal_segments == 1
                and abs(F) > 1e-30):
            ntu = U * math.pi * D_pipe * L / abs(F)
            decay = math.exp(-ntu) if ntu > 0.0 else 1.0
            # A uniform source S shifts the ambient temperature the fluid relaxes
            # to by S/(U·π·D·L); without wall exchange it adds S/(ṁ·cp).
            if ntu > 0.0:
                T_amb_eff = T_amb + S_pipe / (U * math.pi * D_pipe * L)
                rhs_val = (1.0 - decay) * T_amb_eff
            else:
                T_amb_eff = T_amb
                rhs_val = S_pipe / abs(F)
            ntu_paths[pipe_idx] = (ntu, T_amb_eff, S_pipe / abs(F))
            # The sole interior node (offset + 0) represents the pipe exit.
            internal_row = offset
            upstream_col = (node_index[ps.start_node.node_id] if F > 0.0
                            else node_index[ps.end_node.node_id])
            add(internal_row, internal_row, 1.0)
            add(internal_row, upstream_col, -decay)
            rhs[internal_row] = rhs_val
            continue

        # Moukalled source linearisation: S = Sc + Sp·T_P, Sp ≤ 0.
        # The whole wall area U·π·D·L is shared by the n_internal nodes (the
        # junction nodes carry no wall term), so each node gets L/n_internal.
        # Using dx = L/n_segs here would drop 1/n_segs of the pipe area.
        dx_wall = L / n_internal
        Sc = U * math.pi * D_pipe * dx_wall * T_amb
        Sp = -U * math.pi * D_pipe * dx_wall   # negative → strengthens diagonal

        # Fixed power source (heater/cooler): Sc += Q/n_internal, Sp_extra = 0
        if isinstance(ps, HeatSourceState) and n_internal > 0:
            Sc += ps.component.power_w / n_internal
        if S_pipe != 0.0 and n_internal > 0:
            Sc += S_pipe / n_internal

        a_W_conv, a_E_conv = convection_scheme.face_coefficients(F, D)

        junction_start_idx = node_index[ps.start_node.node_id]
        junction_end_idx = node_index[ps.end_node.node_id]

        # Check for Neumann BCs at this pipe's boundary nodes
        start_nid = ps.start_node.node_id
        end_nid = ps.end_node.node_id
        west_neumann = gradient_node_map.get(start_nid)  # None or float
        east_neumann = gradient_node_map.get(end_nid)    # None or float

        for j in range(n_internal):
            row = offset + j

            # ── west neighbour ──
            if j == 0 and west_neumann is not None:
                # Neumann at start node: remove diffusion to start junction
                # a_W_eff = max(F, 0) (convection only), a_P reduced by D
                a_W_j = max(F, 0.0)
                w_col = junction_start_idx
                rhs_extra = -k_fl * A * west_neumann
            elif j == 0:
                a_W_j = a_W_conv
                w_col = junction_start_idx
                rhs_extra = 0.0
            else:
                a_W_j = a_W_conv
                w_col = offset + j - 1
                rhs_extra = 0.0

            # ── east neighbour ──
            if j == n_internal - 1 and east_neumann is not None:
                # Neumann at end node: remove diffusion to end junction
                # a_E_eff = max(-F, 0) (convection only), a_P reduced by D
                a_E_j = max(-F, 0.0)
                e_col = junction_end_idx
                rhs_extra += k_fl * A * east_neumann
            elif j == n_internal - 1:
                a_E_j = a_E_conv
                e_col = junction_end_idx
            else:
                a_E_j = a_E_conv
                e_col = offset + j + 1

            # a_P is always consistent with a_W_j + a_E_j (Neumann removes D from face)
            a_P_j = a_W_j + a_E_j - Sp

            add(row, row, a_P_j)
            rhs[row] += Sc + rhs_extra
            add(row, w_col, -a_W_j)
            add(row, e_col, -a_E_j)

    # ── junction node equations ───────────────────────────────────────────
    inflow: dict[int, list[tuple[int, float]]] = {i: [] for i in range(N_nodes)}
    outflow_total: dict[int, float] = dict.fromkeys(range(N_nodes), 0.0)
    # For compressible fittings the exit temperature T_out ≠ T_in (JT effect).
    # We cannot express T_out as a matrix column, so we put ṁ·Cp·T_out on the
    # RHS directly.  explicit_inflow_total tracks ṁ·Cp so the total_in check
    # does not incorrectly flag the downstream node as degenerate.
    explicit_inflow_total: dict[int, float] = dict.fromkeys(range(N_nodes), 0.0)
    explicit_rhs: dict[int, float] = dict.fromkeys(range(N_nodes), 0.0)

    # Enthalpy-consistent mixing.  Each link transports ṁ·cp_link·T, and cp
    # differs from link to link (pressure, temperature, phase split), so a
    # mixing balance of ṁ·cp·T makes T jump at a node where the fluid merely
    # passes from one pipe into the next.  For fluids with an enthalpy model
    # each stream is written as h ≈ cp·T + [h(T*) − cp·T*] about its previous
    # temperature T*; the bracketed offsets go to the right-hand side, so the
    # junction conserves real enthalpy, Σ ṁ·h_in = Σ ṁ·h_out, at convergence.
    use_offsets = pipe_enthalpy_terms and _has_enthalpy(fluid_model)

    def _offset(link, P: float | None, T_star: float, cp_w: float) -> float:
        if not use_offsets or P is None:
            return 0.0
        h = fluid_enthalpy(fluid_model, link, float(P), T_star)
        return 0.0 if h is None else h - cp_w * T_star

    # Streams entering each node (link, |ṁ|, cp): at a node with no outflow
    # (a sink) the fluid leaving the network is the arriving fluid at the node
    # temperature, so its offset balances the arriving one.
    arriving: dict[int, list[tuple[Any, float, float]]] = {i: [] for i in range(N_nodes)}

    for pipe_idx, ps in enumerate(pipe_states):
        n_segs = max(ps.component.n_thermal_segments, 2)
        n_internal = n_segs - 1
        offset = pipe_internal_offset[pipe_idx]

        T_repr = ps.temperature_c if ps.temperature_c is not None else T_ref
        _ls = _LinkAtTemperature(ps, T_repr)
        cp = fluid_model.specific_heat_for_link(_ls)
        mdot = float(ps.mass_flow_kg_per_s)
        abs_mdot_cp = abs(mdot) * cp

        junction_start_row = node_index[ps.start_node.node_id]
        junction_end_row = node_index[ps.end_node.node_id]

        if abs_mdot_cp < 1e-30:
            continue

        if mdot >= 0.0:
            upstream_col = offset + n_internal - 1 if n_internal > 0 else junction_start_row
            inflow[junction_end_row].append((upstream_col, abs_mdot_cp))
            outflow_total[junction_start_row] += abs_mdot_cp
            up_node, down_node, up_row, down_row = ps.start_node, ps.end_node, junction_start_row, junction_end_row
        else:
            upstream_col = offset + 0 if n_internal > 0 else junction_end_row
            inflow[junction_start_row].append((upstream_col, abs_mdot_cp))
            outflow_total[junction_end_row] += abs_mdot_cp
            up_node, down_node, up_row, down_row = ps.end_node, ps.start_node, junction_end_row, junction_start_row
        if use_offsets:
            T_up_star = up_node.temperature_c if up_node.temperature_c is not None else T_ref
            T_ex = getattr(ps, "exit_temperature_c", None)
            T_ex_star = T_ex if T_ex is not None else T_up_star
            explicit_rhs[up_row] += abs(mdot) * _offset(ps, up_node.pressure_pa, T_up_star, cp)
            explicit_rhs[down_row] -= abs(mdot) * _offset(ps, down_node.pressure_pa, T_ex_star, cp)
            arriving[down_row].append((ps, abs(mdot), cp))

    # ── passthrough devices (fittings, pumps) ────────────────────────────
    # Fittings — adiabatic, no shaft work:
    #   Incompressible: T_out = T_in (implicit matrix term).
    #   Compressible:   isenthalpic, h(T_in,P_in) = h(T_out,P_out).
    # Pumps — shaft work w = ΔP/ρ (J/kg) added to the fluid:
    #   Incompressible: T_out = T_in + ΔP/(ρ·Cp).
    #   Compressible:   h_out = h(T_in,P_in) + ΔP/ρ, then invert for T_out.
    # In all cases where T_out ≠ T_in the downstream contribution is injected
    # as an explicit RHS term so the mixing equation remains valid.
    for dev in network_state.components:
        dt = device_thermal_state(dev, fluid_model, T_ref)
        if dt is None:
            continue
        up_idx = node_index[dt.up_node.node_id]
        down_idx = node_index[dt.down_node.node_id]
        abs_mdot_cp = abs(dt.mdot) * dt.cp_up
        if dt.kind == "passthrough":
            # Incompressible fitting, or a well kept at reservoir temperature:
            # T_out = T_in (implicit matrix term).  With an enthalpy model the
            # stream's offsets are added like for any link; the enthalpy the
            # device changes at constant T is heat exchanged with its
            # surroundings (the formation, for a well) and is reported by the
            # global balance as isothermal-device heat.
            inflow[down_idx].append((up_idx, abs_mdot_cp))
            outflow_total[up_idx] += abs_mdot_cp
            if use_offsets:
                m_abs = abs(dt.mdot)
                explicit_rhs[up_idx] += m_abs * _offset(dev, dt.up_node.pressure_pa, dt.T_up, dt.cp_up)
                explicit_rhs[down_idx] -= m_abs * _offset(dev, dt.down_node.pressure_pa, dt.T_up, dt.cp_up)
                arriving[down_idx].append((dev, m_abs, dt.cp_up))
        else:
            abs_mdot_cp_out = abs(dt.mdot) * dt.cp_out
            outflow_total[up_idx]           += abs_mdot_cp
            explicit_inflow_total[down_idx] += abs_mdot_cp_out
            explicit_rhs[down_idx]          -= abs_mdot_cp_out * dt.T_out
            if use_offsets:
                m_abs = abs(dt.mdot)
                explicit_rhs[up_idx] += m_abs * _offset(dev, dt.up_node.pressure_pa, dt.T_up, dt.cp_up)
                explicit_rhs[down_idx] -= m_abs * _offset(dev, dt.down_node.pressure_pa, dt.T_out, dt.cp_out)
                arriving[down_idx].append((dev, m_abs, dt.cp_out))

    def _sink_offset(nid: int, row: int) -> float:
        """Offset of the fluid leaving the network at a node with no outflow."""
        if not use_offsets or outflow_total[row] > 1e-30:
            return 0.0
        node = network_state.nodes[nid]
        T_star = node.temperature_c if node.temperature_c is not None else T_ref
        return sum(m * _offset(link, node.pressure_pa, T_star, cp_w) for link, m, cp_w in arriving[row])

    thermal_inlet_map = {}
    for node in network_state.nodes.values():
        if node.is_thermal_inlet and node.temperature_c is not None:
            thermal_inlet_map[node.node_id] = node.temperature_c

    for nid in sorted_node_ids:
        row = node_index[nid]

        if nid in thermal_inlet_map:
            # Dirichlet: T = T_prescribed
            add(row, row, 1.0)
            rhs[row] = thermal_inlet_map[nid]

        elif nid in gradient_node_map:
            # Neumann (zero or fixed gradient): extrapolate from first interior node
            # Nodes are spaced dx apart: T_junction = T_adjacent ± g * dx
            pipe_info = gradient_node_pipe.get(nid)
            if pipe_info is not None:
                _pi, is_start, dx_p, offset_p, n_internal_p = pipe_info
                g = gradient_node_map[nid]
                if is_start:
                    # Start junction: T_start = T_j0 - g * dx
                    # (extrapolate backwards: T at x=0 given dT/dx=g and T at x=dx)
                    j0_col = offset_p
                    add(row, row, 1.0)
                    add(row, j0_col, -1.0)
                    rhs[row] = -g * dx_p
                else:
                    # End junction: T_end = T_last + g * dx
                    # (extrapolate forwards: T at x=L given dT/dx=g and T at x=L-dx)
                    last_col = offset_p + n_internal_p - 1
                    add(row, row, 1.0)
                    add(row, last_col, -1.0)
                    rhs[row] = g * dx_p
            else:
                # No interior nodes available: fall back to mixing or 20 °C
                total_in = sum(w for _, w in inflow[row]) + explicit_inflow_total[row]
                total_out = outflow_total[row]
                if total_in < 1e-30:
                    add(row, row, 1.0)
                    rhs[row] = T_ref
                else:
                    denom = total_out if total_out > 1e-30 else total_in
                    add(row, row, -denom)
                    for col, weight in inflow[row]:
                        add(row, col, weight)
                    rhs[row] = explicit_rhs[row] + _sink_offset(nid, row)

        else:
            # Mixing equation (default for junctions and unspecified boundary nodes)
            total_in = sum(w for _, w in inflow[row]) + explicit_inflow_total[row]
            total_out = outflow_total[row]
            if total_in < 1e-30:
                add(row, row, 1.0)
                rhs[row] = T_ref
            else:
                denom = total_out if total_out > 1e-30 else total_in
                add(row, row, -denom)
                for col, weight in inflow[row]:
                    add(row, col, weight)
                rhs[row] = explicit_rhs[row] + _sink_offset(nid, row)

    # ── solve ─────────────────────────────────────────────────────────────
    A_mat = sp.csr_matrix((vals, (rows, cols)), shape=(N_total, N_total))
    T_vec = spla.spsolve(A_mat, rhs)

    # Mean temperature of each pipe's internal FV nodes — more accurate than
    # averaging the two boundary junction temperatures for property evaluation.
    component_mean_temps: dict[int, float] = {}
    for pipe_idx, ps in enumerate(pipe_states):
        n_segs = max(ps.component.n_thermal_segments, 2)
        n_internal = n_segs - 1
        offset = pipe_internal_offset[pipe_idx]
        if n_internal > 0:
            # Temperature of the fluid leaving the pipe: the last interior node
            # in the flow direction (for the NTU path, the sole interior node).
            exit_row = offset + n_internal - 1 if ps.mass_flow_kg_per_s >= 0.0 else offset
            ps.exit_temperature_c = float(T_vec[exit_row])
            if pipe_idx in ntu_paths:
                # Analytical path: the interior node is the exit, so use the exact
                # length average T_amb + (T_in − T_amb)·(1 − e^−NTU)/NTU for
                # property evaluation instead of the (colder) exit temperature.
                ntu, T_amb_p, dT_src = ntu_paths[pipe_idx]
                up_id = ps.start_node.node_id if ps.mass_flow_kg_per_s >= 0.0 else ps.end_node.node_id
                T_in_p = float(T_vec[node_index[up_id]])
                if ntu > 1e-12:
                    factor = (1.0 - math.exp(-ntu)) / ntu
                    component_mean_temps[pipe_idx] = T_amb_p + (T_in_p - T_amb_p) * factor
                else:
                    # No wall exchange: the source makes the profile linear.
                    component_mean_temps[pipe_idx] = T_in_p + 0.5 * dT_src
            else:
                component_mean_temps[pipe_idx] = float(np.mean(T_vec[offset:offset + n_internal]))
        else:
            si = node_index[ps.start_node.node_id]
            ei = node_index[ps.end_node.node_id]
            component_mean_temps[pipe_idx] = float(0.5 * (T_vec[si] + T_vec[ei]))

    node_temps = {nid: float(T_vec[node_index[nid]]) for nid in sorted_node_ids}
    return node_temps, component_mean_temps


@dataclass(frozen=True)
class DeviceThermalState:
    """Inlet/outlet thermal state of a pump or fitting, as used by the solver.

    kind:
      ``"pump"``         shaft work w = ΔP/ρ added (T_out = T_up + w/cp, or
                          h_out = h(T_up, P_up) + w for fluids with enthalpy)
      ``"isenthalpic"``  compressible fitting, h(T_out, P_down) = h(T_up, P_up)
      ``"passthrough"``  incompressible fitting, T_out = T_up
    """

    mdot: float
    up_node: Any
    down_node: Any
    T_up: float
    cp_up: float
    T_out: float
    cp_out: float
    kind: str


def device_thermal_state(dev, fluid_model, T_ref: float) -> DeviceThermalState | None:
    """Thermal state of a pump/fitting link (None for other links or no flow).

    Shared by the energy assembly and the global energy balance so both use
    exactly the same outlet temperature and heat capacities.
    """
    if not isinstance(dev, (FittingState, PumpState)):
        return None
    mdot = float(dev.mass_flow_kg_per_s)
    if abs(mdot) < 1e-30:
        return None
    # Wells (IPR and IPR + tubing) keep the reservoir temperature until the
    # tubing thermal model exists: an isenthalpic expansion without the
    # potential energy of a 2 km lift would be worse than no change.
    is_well = (getattr(dev.component, "is_well_model", False)
               or getattr(dev.component, "is_well_inflow", False))

    up_node   = dev.start_node if mdot >= 0.0 else dev.end_node
    down_node = dev.end_node   if mdot >= 0.0 else dev.start_node
    T_up = up_node.temperature_c if up_node.temperature_c is not None else T_ref
    cp = fluid_model.specific_heat_for_link(_LinkAtTemperature(dev, T_up))

    _p_up   = getattr(up_node,   "pressure_pa", None)
    _p_down = getattr(down_node, "pressure_pa", None)
    _p_ref  = getattr(fluid_model, "reference_pressure_pa", 101_325.0)
    P_up   = _p_up   if _p_up   is not None else _p_ref
    P_down = _p_down if _p_down is not None else _p_ref

    h_up = None if is_well else fluid_enthalpy(fluid_model, dev, P_up, T_up)
    if isinstance(dev, PumpState):
        rho = fluid_model.density_for_link(dev)
        w_shaft = (P_down - P_up) / rho
        if h_up is not None:
            T_out = _temperature_from_enthalpy(fluid_model, dev, h_up + w_shaft, P_down, T_up)
            cp_out = _specific_heat_at(fluid_model, dev, P_down, T_out)
        else:
            T_out = T_up + w_shaft / cp
            cp_out = cp
        return DeviceThermalState(mdot, up_node, down_node, T_up, cp, T_out, cp_out, "pump")
    if h_up is not None:
        T_out = _temperature_from_enthalpy(fluid_model, dev, h_up, P_down, T_up)
        cp_out = _specific_heat_at(fluid_model, dev, P_down, T_out)
        return DeviceThermalState(mdot, up_node, down_node, T_up, cp, T_out, cp_out, "isenthalpic")
    return DeviceThermalState(mdot, up_node, down_node, T_up, cp, T_up, cp, "passthrough")


def _temperature_from_enthalpy(fluid_model, link, h: float, P: float, T_guess: float) -> float:
    if hasattr(fluid_model, "temperature_from_enthalpy_for_link"):
        return fluid_model.temperature_from_enthalpy_for_link(link, h, P, T_guess_c=T_guess)
    return fluid_model.temperature_from_enthalpy(h, P, T_guess_c=T_guess)


def _specific_heat_at(fluid_model, link, P: float, T: float) -> float:
    if hasattr(fluid_model, "specific_heat_at"):
        return fluid_model.specific_heat_at(link, P, T)
    return fluid_model._specific_heat_fn(P, T)


class _LinkAtTemperature:
    """A link state seen at a given temperature.

    Every attribute (composition ``zs``, end-node pressures, component, …) is
    delegated to the real link state; only ``temperature_c`` is overridden.
    Fluid models whose properties depend on pressure or composition — the
    compositional EOS fluid, the black-oil per-inlet fluids — therefore see the
    actual pipe conditions.  (A temperature-only stand-in made the
    compositional fluid fall back to its default composition at 1 atm.)
    """

    __slots__ = ("_link", "temperature_c")

    def __init__(self, link_state, temperature_c: float) -> None:
        self._link = link_state
        self.temperature_c = temperature_c

    def __getattr__(self, name: str):
        return getattr(self._link, name)
