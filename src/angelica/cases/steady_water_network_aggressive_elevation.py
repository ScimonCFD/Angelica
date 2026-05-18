from __future__ import annotations

from dataclasses import replace

from angelica.core.case import NetworkCase
from angelica.core.components import Pipe

from .steady_water_network import build_steady_water_network_case


def build_steady_water_network_aggressive_elevation_case() -> NetworkCase:
    base_case = build_steady_water_network_case()
    height_pattern_m = {
        (1, 7): -2.0,
        (8, 2): -3.0,
        (2, 9): 4.0,
        (10, 3): 3.0,
        (2, 4): 2.0,
        (4, 3): 5.0,
        (5, 4): -3.0,
        (3, 6): 4.0,
    }

    components = []
    for component in base_case.components:
        if isinstance(component, Pipe):
            components.append(
                replace(
                    component,
                    height_change_m=height_pattern_m.get(
                        (component.start_node, component.end_node),
                        0.0,
                    ),
                )
            )
        else:
            components.append(component)

    return replace(
        base_case,
        name="Steady water network with elevation changes",
        components=tuple(components),
    )
