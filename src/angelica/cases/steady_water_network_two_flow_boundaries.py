from __future__ import annotations

from dataclasses import replace

from angelica.core.case import FlowBoundary, NetworkCase, PressureBoundary

from .steady_water_network import build_steady_water_network_case


def build_steady_water_network_two_flow_boundaries_case() -> NetworkCase:
    base_case = build_steady_water_network_case()
    return replace(
        base_case,
        name="Recovered steady water network with two flow boundaries",
        pressure_inlets=(
            PressureBoundary(node_id=5, pressure_pa=201.3 * 1000.0),
        ),
        pressure_outlets=(),
        flow_inlets=(
            FlowBoundary(node_id=1, mass_flow_kg_per_s=1.705309350650952),
        ),
        flow_outlets=(
            FlowBoundary(node_id=6, mass_flow_kg_per_s=2.5530058971552076),
        ),
    )
