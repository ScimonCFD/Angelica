from __future__ import annotations

from dataclasses import replace

from angelica.core.case import FlowBoundary, NetworkCase, PressureBoundary

from .steady_water_network_aggressive_elevation import (
    build_steady_water_network_aggressive_elevation_case,
)


def build_steady_water_network_aggressive_elevation_outlet_flow_case() -> NetworkCase:
    base_case = build_steady_water_network_aggressive_elevation_case()
    return replace(
        base_case,
        name="Steady water network with elevation changes and outlet flow boundary",
        pressure_inlets=(
            PressureBoundary(node_id=1, pressure_pa=251.3 * 1000.0),
            PressureBoundary(node_id=5, pressure_pa=201.3 * 1000.0),
        ),
        pressure_outlets=(),
        flow_outlets=(
            FlowBoundary(node_id=6, mass_flow_kg_per_s=1.7716939852388593),
        ),
        initial_node_pressures_pa={
            1: 251300.0,
            2: 273884.3745906513,
            3: 152253.6637384373,
            4: 227400.6550512336,
            5: 201300.0,
            6: 101300.0,
            7: 251398.75667011723,
            8: 251008.50095061227,
            9: 217825.8550302541,
            10: 215421.51410923176,
        },
    )
