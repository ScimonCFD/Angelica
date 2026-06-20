from __future__ import annotations

from angelica.core.case import NetworkCase, PressureBoundary, ThermalBoundary
from angelica.core.components import Pipe
from angelica.properties.dead_oil import build_thermal_dead_oil

_API_GRAVITY = 32.0   # °API


def build_crude_oil_pipeline_thermal_case() -> NetworkCase:
    """Branched crude oil gathering pipeline with heat loss — 32°API crude.

    Oil enters at 80 °C to benefit from reduced viscosity during transport.
    The pipeline loses heat to a 15 °C ambient through light insulation
    (U = 2 W/(m²·K)), so viscosity rises along the pipeline.  Fluid
    properties are temperature-dependent via:

        Density     : Beggs & Robinson standard-condition value (constant).
        Viscosity   : Beggs & Robinson (1975) dead-oil correlation.
        Specific heat: Watson-Nelson (API Technical Data Book).
        Conductivity : Cragoe (1929).

    Node IDs:
        1 — Pressure inlet      (800 000 Pa, T_in = 80 °C)
        2 — Interior junction
        3 — Pressure outlet A   (101 325 Pa)
        4 — Pressure outlet B   (101 325 Pa)

    Pipes (carbon-steel, ε = 4.6 × 10⁻⁵ m, U = 2 W/(m²·K), T_amb = 15 °C):
        Pipe 1:  1 → 2,  D = 0.10 m, L = 2 000 m  (trunk)
        Pipe 2:  2 → 3,  D = 0.08 m, L = 1 000 m  (branch A)
        Pipe 3:  2 → 4,  D = 0.06 m, L =   600 m  (branch B)
    """
    fluid = build_thermal_dead_oil(_API_GRAVITY)
    roughness = 4.6e-5   # carbon-steel, m
    U = 2.0              # W/(m²·K) — light insulation
    T_amb = 15.0         # °C
    n_seg = 30           # FV segments per pipe

    return NetworkCase(
        name="Crude Oil Pipeline — thermal (32°API, T_in=80°C)",
        fluid_model=fluid,
        pressure_inlets=(
            PressureBoundary(node_id=1, pressure_pa=800_000.0),
        ),
        pressure_outlets=(
            PressureBoundary(node_id=3, pressure_pa=101_325.0),
            PressureBoundary(node_id=4, pressure_pa=101_325.0),
        ),
        thermal_inlets=(
            ThermalBoundary(node_id=1, temperature_c=80.0),
        ),
        components=(
            Pipe(1, 2, diameter_m=0.10, length_m=2000.0,
                 absolute_roughness_m=roughness,
                 heat_transfer_coefficient_w_per_m2k=U,
                 ambient_temperature_c=T_amb,
                 n_thermal_segments=n_seg,
                 component_id="trunk"),
            Pipe(2, 3, diameter_m=0.08, length_m=1000.0,
                 absolute_roughness_m=roughness,
                 heat_transfer_coefficient_w_per_m2k=U,
                 ambient_temperature_c=T_amb,
                 n_thermal_segments=n_seg,
                 component_id="branch_a"),
            Pipe(2, 4, diameter_m=0.06, length_m=600.0,
                 absolute_roughness_m=roughness,
                 heat_transfer_coefficient_w_per_m2k=U,
                 ambient_temperature_c=T_amb,
                 n_thermal_segments=n_seg,
                 component_id="branch_b"),
        ),
        node_ids=(1, 2, 3, 4),
    )
