from __future__ import annotations

from angelica.core.case import NetworkCase, PressureBoundary
from angelica.core.components import Pipe
from angelica.properties.single_component import SingleComponentFluid


def build_laminar_parallel_pipes_case() -> NetworkCase:
    """Three parallel pipes between pressure inlet and outlet — laminar regime.

    Oil fluid (ρ=880 kg/m³, μ=0.1 Pa·s) flows through three smooth parallel
    pipes from node 1 (P_in=50 000 Pa) to node 2 (P_out=0 Pa).  All pipes are
    smooth (roughness=0) and all operate well within the laminar regime (Re<180).

    Node IDs:
        1 — Pressure inlet  (50 000 Pa)
        2 — Pressure outlet (0 Pa)

    Pipes (all from node 1 to node 2):
        Pipe 1:  D=0.020 m, L=10 m
        Pipe 2:  D=0.015 m, L= 8 m
        Pipe 3:  D=0.025 m, L=12 m

    Analytical Poiseuille solution:  Q = π D⁴ ΔP / (128 μ L)
    """
    return NetworkCase(
        name="Laminar Poiseuille Benchmark",
        fluid_model=SingleComponentFluid(
            density_kg_per_m3=880.0,
            viscosity_pa_s=0.1,
        ),
        pressure_inlets=(
            PressureBoundary(node_id=1, pressure_pa=50000.0),
        ),
        pressure_outlets=(
            PressureBoundary(node_id=2, pressure_pa=0.0),
        ),
        components=(
            Pipe(1, 2, diameter_m=0.020, length_m=10.0, absolute_roughness_m=0.0,
                 height_change_m=0.0, component_id="pipe_1"),
            Pipe(1, 2, diameter_m=0.015, length_m=8.0,  absolute_roughness_m=0.0,
                 height_change_m=0.0, component_id="pipe_2"),
            Pipe(1, 2, diameter_m=0.025, length_m=12.0, absolute_roughness_m=0.0,
                 height_change_m=0.0, component_id="pipe_3"),
        ),
        node_ids=(1, 2),
    )
