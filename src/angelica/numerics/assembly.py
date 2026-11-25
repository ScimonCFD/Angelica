from __future__ import annotations

import numpy as np
from scipy.sparse import lil_matrix


def assemble_pressure_system(network_state, couplings: list[float]):
    node_ids = sorted(network_state.nodes)
    node_index = {node_id: idx for idx, node_id in enumerate(node_ids)}
    n = len(node_ids)
    A = lil_matrix((n, n), dtype=float)
    b = np.zeros(n, dtype=float)

    separators = getattr(network_state, "separators", None) or {}
    # A controlling regulator limited by its wide-open capacity in this
    # iteration acts as the open valve (K-factor coupling, flow as computed).
    regulators = [r for r in (getattr(network_state, "regulators", None) or ())
                  if r.regulator_state == "controlling" and not r.regulator_capped]
    for link_state, coupling in zip(network_state.components, couplings):
        link_state.coupling = coupling
    if separators:
        _separator_outlet_flows(separators)
    if regulators:
        _regulator_flows(regulators)

    for link_state, coupling in zip(network_state.components, couplings):
        i = node_index[link_state.start_node.node_id]
        j = node_index[link_state.end_node.node_id]

        b[i] -= link_state.mass_flow_kg_per_s
        b[j] += link_state.mass_flow_kg_per_s
        A[i, j] += coupling
        A[j, i] += coupling
        A[i, i] -= coupling
        A[j, j] -= coupling

    if separators:
        _separator_outlet_jacobian(separators, node_index, A)
    if regulators:
        _regulator_jacobian(regulators, node_index, A, separators)

    for node_id in node_ids:
        node = network_state.nodes[node_id]
        if node.is_pressure_boundary or node.regulated:
            idx = node_index[node_id]
            A.rows[idx] = []
            A.data[idx] = []
            A[idx, idx] = 1.0
            b[idx] = 0.0
        else:
            idx = node_index[node_id]
            if node.prescribed_mass_flow_kg_per_s is not None:
                if node.is_inlet:
                    b[idx] += node.prescribed_mass_flow_kg_per_s
                else:
                    b[idx] -= node.prescribed_mass_flow_kg_per_s
            b[idx] += node.source_mass_flow_kg_per_s      # internal source (lift gas)

    return A.tocsr(), b


# ── separator outlets ─────────────────────────────────────────────────────────
#
# A ``SeparatorOutlet`` carries the fraction f of its separator's feed F (the
# flows arriving through the node's other links): ṁ_out = f·F.  Treating it
# as a fixed flow from the previous iteration makes the pressure correction
# ignore that the outlet follows the feed; with f close to 1 (an LP separator
# whose feed is nearly all liquid) that diverges.  So the outlet flow is set
# from the current feed, and its dependence on the feed links enters the
# matrix: δṁ_out = f·Σ_k |c_k|·(δP_k − δP_S) over the links k flowing into
# the separator S.  Row n of the system reads −Σ δ(net inflow into n) = net
# inflow into n.

def _feed_links_in(outlet):
    node = outlet.start_node
    for link in outlet.separator_feed_links:
        q = float(link.mass_flow_kg_per_s)
        if link.end_node is node and q > 0.0:
            yield link, link.start_node
        elif link.start_node is node and q < 0.0:
            yield link, link.end_node


def _separator_outlet_flows(separators) -> None:
    from angelica.closures.separator import separator_outlet_mass_flow

    for outlets in separators.values():
        for o in outlets:
            o.mass_flow_kg_per_s = separator_outlet_mass_flow(o)


def _separator_outlet_jacobian(separators, node_index, A) -> None:
    for outlets in separators.values():
        for o in outlets:
            f = o.separator_fraction
            if f <= 0.0:
                continue
            s = node_index[o.start_node.node_id]
            out = node_index[o.end_node.node_id]
            for link, other in _feed_links_in(o):
                g = f * abs(float(link.coupling))
                k = node_index[other.node_id]
                A[s, k] += g
                A[s, s] -= g
                A[out, k] -= g
                A[out, s] += g


# ── pressure regulators ───────────────────────────────────────────────────────
#
# A controlling ``PressureRegulator`` holds its controlled node C at the set
# point (a fixed-pressure row) and carries the flow that balances C's other
# links: ṁ_r = Σ_k (flow from C through k) for a reducing regulator (C its
# end node), the supply arriving at C for a back-pressure regulator (C its
# start node).  With δP_C = 0 that flow moves with the pressures of C's
# neighbours: δṁ_r = ∓Σ_k |c_k|·δP_k, which enters the row of the
# regulator's other node O (net inflow into O changes by +|c_k|·δP_k in
# both modes).  Lagging ṁ_r instead is the separator-outlet mistake of
# v1.26.0.

def _regulator_flows(regulators) -> None:
    from angelica.closures.regulator import regulator_mass_flow

    for r in regulators:
        r.mass_flow_kg_per_s = regulator_mass_flow(r)


def _regulator_jacobian(regulators, node_index, A, separators=None) -> None:
    from angelica.closures.regulator import controlled_node

    for r in regulators:
        c = controlled_node(r)
        other = r.start_node if c is r.end_node else r.end_node
        o = node_index[other.node_id]
        # At a separator node the liquid outlets take the fraction F of
        # what arrives, so only (1 − F) of a change in a feed link reaches
        # the regulator (a back-pressure regulator on a separator's gas line).
        outlets = (separators or {}).get(c.node_id, ())
        scale = max(1.0 - sum(o_.separator_fraction for o_ in outlets), 0.0)
        for link in r.regulator_links:
            if any(link is o_ for o_ in outlets):
                continue
            k_node = link.start_node if link.end_node is c else link.end_node
            A[o, node_index[k_node.node_id]] -= scale * abs(float(link.coupling))
