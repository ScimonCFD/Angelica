from __future__ import annotations

import math
import warnings

import numpy as np

from angelica.closures.convection_scheme import ConvectionScheme
from angelica.core.state import HeatSourceState, NetworkState, PipeState


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


def solve_energy_system(
    network_state: NetworkState,
    fluid_model,
    convection_scheme: ConvectionScheme,
    T_ref: float = 20.0,
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
            U = pipe.heat_transfer_coefficient_w_per_m2k
            T_amb = pipe.ambient_temperature_c

        T_repr = ps.temperature_c if hasattr(ps, "temperature_c") and ps.temperature_c is not None else T_ref
        _ls = _TempCarrier(T_repr)

        cp = fluid_model.specific_heat_for_link(_ls)
        k_fl = fluid_model.thermal_conductivity_for_link(_ls)
        mdot = float(ps.mass_flow_kg_per_s)

        F = mdot * cp          # convective strength (W/K), signed
        D = k_fl * A / dx      # diffusive conductance (W/K)

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
            # The sole interior node (offset + 0) represents the pipe exit.
            internal_row = offset
            upstream_col = (node_index[ps.start_node.node_id] if F > 0.0
                            else node_index[ps.end_node.node_id])
            add(internal_row, internal_row, 1.0)
            add(internal_row, upstream_col, -decay)
            rhs[internal_row] = (1.0 - decay) * T_amb
            continue

        # Moukalled source linearisation: S = Sc + Sp·T_P, Sp ≤ 0
        Sc = U * math.pi * D_pipe * dx * T_amb
        Sp = -U * math.pi * D_pipe * dx   # negative → strengthens diagonal

        # Fixed power source (heater/cooler): Sc += Q/n_internal, Sp_extra = 0
        n_internal = n_segs - 1
        if isinstance(ps, HeatSourceState) and n_internal > 0:
            Sc += ps.component.power_w / n_internal

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
    outflow_total: dict[int, float] = {i: 0.0 for i in range(N_nodes)}

    for pipe_idx, ps in enumerate(pipe_states):
        n_segs = max(ps.component.n_thermal_segments, 2)
        n_internal = n_segs - 1
        offset = pipe_internal_offset[pipe_idx]

        T_repr = ps.temperature_c if hasattr(ps, "temperature_c") and ps.temperature_c is not None else T_ref
        _ls = _TempCarrier(T_repr)
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
        else:
            upstream_col = offset + 0 if n_internal > 0 else junction_end_row
            inflow[junction_start_row].append((upstream_col, abs_mdot_cp))
            outflow_total[junction_end_row] += abs_mdot_cp

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
                total_in = sum(w for _, w in inflow[row])
                total_out = outflow_total[row]
                if total_in < 1e-30:
                    add(row, row, 1.0)
                    rhs[row] = T_ref
                else:
                    denom = total_out if total_out > 1e-30 else total_in
                    add(row, row, -denom)
                    for col, weight in inflow[row]:
                        add(row, col, weight)
                    rhs[row] = 0.0

        else:
            # Mixing equation (default for junctions and unspecified boundary nodes)
            total_in = sum(w for _, w in inflow[row])
            total_out = outflow_total[row]
            if total_in < 1e-30:
                add(row, row, 1.0)
                rhs[row] = T_ref
            else:
                denom = total_out if total_out > 1e-30 else total_in
                add(row, row, -denom)
                for col, weight in inflow[row]:
                    add(row, col, weight)
                rhs[row] = 0.0

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
            component_mean_temps[pipe_idx] = float(np.mean(T_vec[offset:offset + n_internal]))
        else:
            si = node_index[ps.start_node.node_id]
            ei = node_index[ps.end_node.node_id]
            component_mean_temps[pipe_idx] = float(0.5 * (T_vec[si] + T_vec[ei]))

    node_temps = {nid: float(T_vec[node_index[nid]]) for nid in sorted_node_ids}
    return node_temps, component_mean_temps


class _TempCarrier:
    """Minimal stand-in for link_state when only temperature_c is needed."""

    __slots__ = ("temperature_c",)

    def __init__(self, temperature_c: float) -> None:
        self.temperature_c = temperature_c
