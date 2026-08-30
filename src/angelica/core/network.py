from __future__ import annotations

from .case import NetworkCase
from .components import Fitting, HeatSource, Pipe, Pump
from .state import FittingState, HeatSourceState, NetworkState, NodeState, PipeState, PressureChangerState, PumpState


def build_network_state(case: NetworkCase) -> NetworkState:
    nodes = {node_id: NodeState(node_id=node_id) for node_id in case.all_node_ids()}

    for node_id, pressure_pa in case.initial_node_pressures_pa.items():
        if node_id in nodes:
            nodes[node_id].pressure_pa = pressure_pa

    for boundary in case.pressure_inlets:
        node = nodes[boundary.node_id]
        node.pressure_pa = boundary.pressure_pa
        node.is_boundary = True
        node.is_inlet = True
        node.is_pressure_boundary = True

    for boundary in case.pressure_outlets:
        node = nodes[boundary.node_id]
        node.pressure_pa = boundary.pressure_pa
        node.is_boundary = True
        node.is_inlet = False
        node.is_pressure_boundary = True

    for fb in case.flow_inlets:
        node = nodes[fb.node_id]
        node.is_boundary = True
        node.is_inlet = True
        node.prescribed_mass_flow_kg_per_s = fb.mass_flow_kg_per_s

    for fb in case.flow_outlets:
        node = nodes[fb.node_id]
        node.is_boundary = True
        node.is_inlet = False
        node.prescribed_mass_flow_kg_per_s = fb.mass_flow_kg_per_s

    components: list[PressureChangerState] = []
    for component in case.components:
        start_node = nodes[component.start_node]
        end_node = nodes[component.end_node]
        if isinstance(component, Pipe):
            components.append(PipeState(component=component, start_node=start_node, end_node=end_node))
        elif isinstance(component, Fitting):
            components.append(FittingState(component=component, start_node=start_node, end_node=end_node))
        elif isinstance(component, Pump):
            components.append(PumpState(component=component, start_node=start_node, end_node=end_node))
        elif isinstance(component, HeatSource):
            components.append(HeatSourceState(component=component, start_node=start_node, end_node=end_node))
        else:
            raise TypeError(f"Unsupported component type: {type(component).__name__}")

    _check_network_connectivity(nodes, components)
    return NetworkState(nodes=nodes, components=components)


def _check_network_connectivity(nodes: dict, components: list) -> None:
    if not nodes:
        return
    adjacency: dict[int, set[int]] = {nid: set() for nid in nodes}
    for comp in components:
        a, b = comp.start_node.node_id, comp.end_node.node_id
        adjacency[a].add(b)
        adjacency[b].add(a)
    start = next(iter(nodes))
    visited: set[int] = set()
    stack = [start]
    while stack:
        nid = stack.pop()
        if nid in visited:
            continue
        visited.add(nid)
        stack.extend(adjacency[nid] - visited)
    unreachable = set(nodes) - visited
    if unreachable:
        raise ValueError(
            f"Network is not fully connected. Unreachable nodes: {sorted(unreachable)}"
        )
