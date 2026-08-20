from __future__ import annotations

import numpy as np
from scipy.sparse import lil_matrix


def assemble_pressure_system(network_state, couplings: list[float]):
    node_ids = sorted(network_state.nodes)
    node_index = {node_id: idx for idx, node_id in enumerate(node_ids)}
    n = len(node_ids)
    A = lil_matrix((n, n), dtype=float)
    b = np.zeros(n, dtype=float)

    for link_state, coupling in zip(network_state.components, couplings):
        i = node_index[link_state.start_node.node_id]
        j = node_index[link_state.end_node.node_id]
        link_state.coupling = coupling

        b[i] -= link_state.mass_flow_kg_per_s
        b[j] += link_state.mass_flow_kg_per_s
        A[i, j] += coupling
        A[j, i] += coupling
        A[i, i] -= coupling
        A[j, j] -= coupling

    for node_id in node_ids:
        node = network_state.nodes[node_id]
        if node.is_pressure_boundary:
            idx = node_index[node_id]
            A.rows[idx] = []
            A.data[idx] = []
            A[idx, idx] = 1.0
            b[idx] = 0.0
        elif node.prescribed_mass_flow_kg_per_s is not None:
            idx = node_index[node_id]
            if node.is_inlet:
                b[idx] += node.prescribed_mass_flow_kg_per_s
            else:
                b[idx] -= node.prescribed_mass_flow_kg_per_s

    return A.tocsr(), b
