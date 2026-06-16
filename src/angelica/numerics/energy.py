from __future__ import annotations

import math

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from angelica.closures.convection_scheme import ConvectionScheme
from angelica.core.state import NetworkState, PipeState


def solve_energy_system(
    network_state: NetworkState,
    fluid_model,
    convection_scheme: ConvectionScheme,
) -> dict[int, float]:
    """Solve the steady-state energy equation for the pipe network.

    Returns a dict mapping node_id → temperature_c for all network nodes.

    Algorithm
    ---------
    Global unknowns:
      - T at each network node (N_nodes)
      - T at each internal pipe node (N_i - 1 per pipe with N_i segments)

    Equations:
      - Thermal inlet nodes  → Dirichlet: T = T_prescribed
      - Non-inlet network nodes → mixing: T_J = Σ(ṁ_k · T_k^upstream) / Σṁ_k
      - Pipe internal nodes  → FV convection-diffusion (node-centred, upwind default)
    """
    sorted_node_ids = sorted(network_state.nodes.keys())
    N_nodes = len(sorted_node_ids)
    node_index = {nid: i for i, nid in enumerate(sorted_node_ids)}

    # ── index internal pipe nodes ──────────────────────────────────────────
    # Only PipeState objects get internal thermal nodes.
    pipe_states = [ls for ls in network_state.components if isinstance(ls, PipeState)]

    pipe_internal_offset: list[int] = []
    total_internal = 0
    for ps in pipe_states:
        pipe_internal_offset.append(N_nodes + total_internal)
        n_segs = max(ps.component.n_thermal_segments, 2)
        total_internal += n_segs - 1

    N_total = N_nodes + total_internal

    # ── build sparse matrix and RHS ───────────────────────────────────────
    rows: list[int] = []
    cols: list[int] = []
    vals: list[float] = []
    rhs = np.zeros(N_total)

    def add(r: int, c: int, v: float) -> None:
        rows.append(r)
        cols.append(c)
        vals.append(v)

    # ── FV equations for pipe internal nodes ──────────────────────────────
    for pipe_idx, ps in enumerate(pipe_states):
        n_segs = max(ps.component.n_thermal_segments, 2)
        n_internal = n_segs - 1            # nodes strictly inside the pipe
        offset = pipe_internal_offset[pipe_idx]

        pipe = ps.component
        L = pipe.length_m
        dx = L / n_segs
        A = ps.area_m2
        D_pipe = pipe.diameter_m
        U = pipe.heat_transfer_coefficient_w_per_m2k
        T_amb = pipe.ambient_temperature_c

        # representative temperature for property evaluation
        T_repr = ps.temperature_c if hasattr(ps, "temperature_c") and ps.temperature_c is not None else 20.0
        # Use a dummy link_state-like object that carries temperature
        _ls = _TempCarrier(T_repr)

        cp = fluid_model.specific_heat_for_link(_ls)
        k_fl = fluid_model.thermal_conductivity_for_link(_ls)
        mdot = float(ps.mass_flow_kg_per_s)  # positive = start→end

        F = mdot * cp          # convective strength (W/K), signed
        D = k_fl * A / dx      # diffusive conductance (W/K)

        # source linearisation: S = U · π · D_pipe · dx · (T_amb - T)
        #   → Sc = U πD dx T_amb  (constant part)
        #   → Sp = -U πD dx       (coefficient of T_P, negative → adds to a_P)
        Sc = U * math.pi * D_pipe * dx * T_amb
        Sp = -U * math.pi * D_pipe * dx   # negative

        a_W_conv, a_E_conv = convection_scheme.face_coefficients(F, D)

        # a_P = a_W + a_E - Sp (Sp is negative, so -Sp is positive)
        a_P = a_W_conv + a_E_conv - Sp

        # index helpers
        junction_start_idx = node_index[ps.start_node.node_id]
        junction_end_idx = node_index[ps.end_node.node_id]

        for j in range(n_internal):  # j = 0 .. n_internal-1 (global internal node j+1)
            row = offset + j
            add(row, row, a_P)
            rhs[row] += Sc

            # west neighbour
            if j == 0:
                w_col = junction_start_idx
            else:
                w_col = offset + j - 1
            add(row, w_col, -a_W_conv)

            # east neighbour
            if j == n_internal - 1:
                e_col = junction_end_idx
            else:
                e_col = offset + j + 1
            add(row, e_col, -a_E_conv)

    # ── junction node equations ───────────────────────────────────────────
    # accumulate per-node: sum of (|ṁ| · cp) entering, and which upstream T
    # For each node: inflow_contributions = list of (global_col, weight)
    #                outflow_total = sum of |ṁ| · cp leaving

    inflow: dict[int, list[tuple[int, float]]] = {i: [] for i in range(N_nodes)}
    outflow_total: dict[int, float] = {i: 0.0 for i in range(N_nodes)}

    for pipe_idx, ps in enumerate(pipe_states):
        n_segs = max(ps.component.n_thermal_segments, 2)
        n_internal = n_segs - 1
        offset = pipe_internal_offset[pipe_idx]

        T_repr = ps.temperature_c if hasattr(ps, "temperature_c") and ps.temperature_c is not None else 20.0
        _ls = _TempCarrier(T_repr)
        cp = fluid_model.specific_heat_for_link(_ls)
        mdot = float(ps.mass_flow_kg_per_s)
        abs_mdot_cp = abs(mdot) * cp

        junction_start_row = node_index[ps.start_node.node_id]
        junction_end_row = node_index[ps.end_node.node_id]

        if abs_mdot_cp < 1e-30:
            continue  # no flow — pipe contributes nothing

        if mdot >= 0.0:
            # flow start→end: enters end junction from last internal node
            upstream_col = offset + n_internal - 1 if n_internal > 0 else junction_start_row
            inflow[junction_end_row].append((upstream_col, abs_mdot_cp))
            outflow_total[junction_start_row] += abs_mdot_cp
        else:
            # flow end→start: enters start junction from first internal node
            upstream_col = offset + 0 if n_internal > 0 else junction_end_row
            inflow[junction_start_row].append((upstream_col, abs_mdot_cp))
            outflow_total[junction_end_row] += abs_mdot_cp

    # build junction rows
    thermal_inlet_map = {}
    for node in network_state.nodes.values():
        if node.is_thermal_inlet and node.temperature_c is not None:
            thermal_inlet_map[node.node_id] = node.temperature_c

    for nid in sorted_node_ids:
        row = node_index[nid]
        if nid in thermal_inlet_map:
            # Dirichlet
            add(row, row, 1.0)
            rhs[row] = thermal_inlet_map[nid]
        else:
            total_in = sum(w for _, w in inflow[row])
            total_out = outflow_total[row]
            if total_in < 1e-30:
                # no inflow information — set T = T_amb (fallback)
                add(row, row, 1.0)
                rhs[row] = 20.0
            else:
                # mixing: T_J · total_out = Σ weight · T_upstream
                # equivalently: -total_out · T_J + Σ weight · T_upstream = 0
                # but if total_out == 0 (pure sink): T_J = mixing of inflows
                # → T_J · total_in = Σ weight · T_upstream
                denom = total_out if total_out > 1e-30 else total_in
                add(row, row, -denom)
                for col, weight in inflow[row]:
                    add(row, col, weight)
                rhs[row] = 0.0

    # ── solve ─────────────────────────────────────────────────────────────
    A_mat = sp.csr_matrix((vals, (rows, cols)), shape=(N_total, N_total))
    T_vec = spla.spsolve(A_mat, rhs)

    return {nid: float(T_vec[node_index[nid]]) for nid in sorted_node_ids}


class _TempCarrier:
    """Minimal stand-in for link_state when only temperature_c is needed."""

    __slots__ = ("temperature_c",)

    def __init__(self, temperature_c: float) -> None:
        self.temperature_c = temperature_c
