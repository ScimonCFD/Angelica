from __future__ import annotations

import numpy as np


def assemble_pressure_system(network_state, couplings: list[float]) -> tuple[np.ndarray, np.ndarray]:
    node_ids = sorted(network_state.nodes)
    node_index = {node_id: idx for idx, node_id in enumerate(node_ids)}
    matrix = np.zeros((len(node_ids), len(node_ids)), dtype=float)
    vector = np.zeros(len(node_ids), dtype=float)

    for link_state, coupling in zip(network_state.components, couplings):
        start = node_index[link_state.start_node.node_id]
        end = node_index[link_state.end_node.node_id]
        link_state.coupling = coupling

        vector[start] -= link_state.mass_flow_kg_per_s
        vector[end] += link_state.mass_flow_kg_per_s
        matrix[start, end] += coupling
        matrix[end, start] += coupling
        matrix[start, start] -= coupling
        matrix[end, end] -= coupling

    for node_id in node_ids:
        node = network_state.nodes[node_id]
        if node.is_pressure_boundary:
            idx = node_index[node_id]
            matrix[idx, :] = 0.0
            matrix[idx, idx] = 1.0
            vector[idx] = 0.0
        elif node.prescribed_mass_flow_kg_per_s is not None:
            idx = node_index[node_id]
            if node.is_inlet:
                vector[idx] += node.prescribed_mass_flow_kg_per_s
            else:
                vector[idx] -= node.prescribed_mass_flow_kg_per_s

    return matrix, vector
