from __future__ import annotations

from dataclasses import replace

from angelica.core.case import FlowBoundary, NetworkCase, PressureBoundary

from .steady_water_network import build_steady_water_network_case


def build_steady_water_network_inlet_flow_boundary_case() -> NetworkCase:
    base_case = build_steady_water_network_case()
    return replace(
        base_case,
        name="Recovered steady water network with inlet flow boundary",
        pressure_inlets=(
            PressureBoundary(node_id=1, pressure_pa=251.3 * 1000.0),
        ),
        pressure_outlets=(
            PressureBoundary(node_id=6, pressure_pa=101.3 * 1000.0),
        ),
        flow_inlets=(
            FlowBoundary(node_id=5, mass_flow_kg_per_s=0.8476957104642427),
        ),
        initial_node_pressures_pa={
            1: 251300.0,
            2: 214100.8325677869,
            3: 124591.49451116432,
            4: 185985.3117864101,
            5: 201300.0,
            6: 101300.0,
            7: 223825.66435200584,
            8: 223258.94445045164,
            9: 185651.9362986657,
            10: 181489.28704940667,
        },
    )
