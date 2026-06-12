from __future__ import annotations

from angelica.core.case import NetworkCase, PressureBoundary
from angelica.core.components import Pipe
from angelica.properties import dead_oil_density_kg_per_m3, dead_oil_viscosity_pa_s
from angelica.properties.single_component import SingleComponentFluid

_API_GRAVITY = 32.0   # °API
_TEMPERATURE_C = 65.0  # °C — used only to evaluate viscosity


def build_crude_oil_pipeline_case() -> NetworkCase:
    """Branched crude oil gathering pipeline — 32°API crude at 65 °C.

    Fluid properties are computed from API gravity and temperature via the
    Beggs & Robinson (1975) dead-oil correlation:
        density  ≈ 864.6 kg/m³
        viscosity ≈ 4.25 cP (0.004254 Pa·s)

    Node IDs:
        1 — Pressure inlet      (500 000 Pa)
        2 — Interior junction
        3 — Pressure outlet A   (101 325 Pa)
        4 — Pressure outlet B   (101 325 Pa)

    Pipes (carbon-steel, ε = 4.6 × 10⁻⁵ m):
        Pipe 1:  1 → 2,  D=0.10 m, L=1 000 m
        Pipe 2:  2 → 3,  D=0.08 m, L=500 m
        Pipe 3:  2 → 4,  D=0.06 m, L=300 m
    """
    rho = dead_oil_density_kg_per_m3(_API_GRAVITY)
    mu = dead_oil_viscosity_pa_s(_API_GRAVITY, _TEMPERATURE_C)
    roughness = 4.6e-5  # carbon-steel commercial pipe, m

    return NetworkCase(
        name="Crude Oil Pipeline (32°API, 65°C)",
        fluid_model=SingleComponentFluid(
            density_kg_per_m3=rho,
            viscosity_pa_s=mu,
        ),
        pressure_inlets=(
            PressureBoundary(node_id=1, pressure_pa=500_000.0),
        ),
        pressure_outlets=(
            PressureBoundary(node_id=3, pressure_pa=101_325.0),
            PressureBoundary(node_id=4, pressure_pa=101_325.0),
        ),
        components=(
            Pipe(1, 2, diameter_m=0.10, length_m=1000.0,
                 absolute_roughness_m=roughness, height_change_m=0.0,
                 component_id="pipe_1"),
            Pipe(2, 3, diameter_m=0.08, length_m=500.0,
                 absolute_roughness_m=roughness, height_change_m=0.0,
                 component_id="pipe_2"),
            Pipe(2, 4, diameter_m=0.06, length_m=300.0,
                 absolute_roughness_m=roughness, height_change_m=0.0,
                 component_id="pipe_3"),
        ),
        node_ids=(1, 2, 3, 4),
    )
