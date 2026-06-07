from __future__ import annotations

from angelica.core.case import NetworkCase, PressureBoundary
from angelica.core.components import Pipe
from angelica.properties.single_component import SingleComponentFluid

_RHO = 998.25
_G = 9.81


def build_three_reservoir_junction_case() -> NetworkCase:
    """Classical three-reservoir junction problem.

    Reservoir A (30 m) feeds junction J; J splits to reservoir B (20 m)
    and reservoir C (5 m).  Pipe sizes chosen so that all three
    pipes operate in the turbulent rough regime.

    Node IDs:
        1 — Reservoir A  (source, 30 m head)
        2 — Junction J   (unknown head)
        3 — Reservoir B  (sink,   20 m head)
        4 — Reservoir C  (sink,    5 m head)

    Pipes:
        A→J  D = 0.200 m, L = 1000 m, k = 0.0001 m
        J→B  D = 0.150 m, L =  800 m, k = 0.0001 m
        J→C  D = 0.100 m, L =  600 m, k = 0.0001 m
    """
    return NetworkCase(
        name="Three-reservoir junction",
        fluid_model=SingleComponentFluid(
            density_kg_per_m3=_RHO,
            viscosity_pa_s=0.001,
        ),
        pressure_inlets=(
            PressureBoundary(node_id=1, pressure_pa=_RHO * _G * 30.0),
        ),
        pressure_outlets=(
            PressureBoundary(node_id=3, pressure_pa=_RHO * _G * 20.0),
            PressureBoundary(node_id=4, pressure_pa=_RHO * _G * 5.0),
        ),
        components=(
            Pipe(1, 2, diameter_m=0.200, length_m=1000.0, absolute_roughness_m=0.0001, height_change_m=0.0),
            Pipe(2, 3, diameter_m=0.150, length_m=800.0,  absolute_roughness_m=0.0001, height_change_m=0.0),
            Pipe(2, 4, diameter_m=0.100, length_m=600.0,  absolute_roughness_m=0.0001, height_change_m=0.0),
        ),
        node_ids=(1, 2, 3, 4),
    )
