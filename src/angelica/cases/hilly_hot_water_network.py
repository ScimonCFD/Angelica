from __future__ import annotations

from angelica.core.case import FlowBoundary, NetworkCase, PressureBoundary, ThermalBoundary
from angelica.core.components import Pipe
from angelica.properties.water_liquid import build_water_thermal_fluid


def build_hilly_hot_water_network_case() -> NetworkCase:
    """Hilly hot-water district heating — mixed elevation directions + outflow BC.

    A supply main climbs from a valve station (node 1) to a hilltop junction
    (node 2), then splits into two downhill branches.  Branch A (2→3) descends
    steeply to a fixed-flow extraction point; branch B (2→4) descends gently
    to a pressure-controlled return header.

    All pipe lengths are measured along the pipe axis (actual inclined length),
    so length_m > |height_change_m| in every pipe.

    Elevation convention (relative to node 1):
        Node 1  — base   (supply valve station)
        Node 2  — +40 m  (hilltop junction, ascending trunk)
        Node 3  — −20 m  (consumer A, 60 m below hilltop)
        Node 4  — +10 m  (consumer B, 30 m below hilltop)

    Pipes:
        1 → 2  D=0.15 m  L=150 m  Δz=+40 m  (trunk, ascending)
        2 → 3  D=0.10 m  L=220 m  Δz=−60 m  (branch A, steeply descending)
        2 → 4  D=0.12 m  L=120 m  Δz=−30 m  (branch B, gently descending)

    Boundary conditions:
        Node 1 — pressure inlet  (600 000 Pa, T_in = 85 °C)
        Node 3 — outflow BC      (ṁ_out = 3.0 kg/s, pressure free)
        Node 4 — pressure outlet (200 000 Pa)
    """
    roughness = 4.6e-5   # carbon-steel commercial pipe, m
    U     = 8.0          # W/(m²·K) — moderate insulation
    T_amb = 10.0         # °C
    n_seg = 20

    return NetworkCase(
        name="Hilly hot-water district heating — mixed elevation + outflow BC",
        fluid_model=build_water_thermal_fluid(),
        pressure_inlets=(
            PressureBoundary(node_id=1, pressure_pa=600_000.0),
        ),
        pressure_outlets=(
            PressureBoundary(node_id=4, pressure_pa=200_000.0),
        ),
        flow_outlets=(
            FlowBoundary(node_id=3, mass_flow_kg_per_s=3.0),
        ),
        thermal_inlets=(
            ThermalBoundary(node_id=1, temperature_c=85.0, bc_type="fixed_temperature"),
            ThermalBoundary(node_id=3, bc_type="zero_gradient"),
            ThermalBoundary(node_id=4, bc_type="zero_gradient"),
        ),
        components=(
            Pipe(1, 2,
                 diameter_m=0.15, length_m=150.0,
                 absolute_roughness_m=roughness,
                 height_change_m=+40.0,
                 heat_transfer_coefficient_w_per_m2k=U,
                 ambient_temperature_c=T_amb,
                 n_thermal_segments=n_seg,
                 component_id="trunk_ascending"),
            Pipe(2, 3,
                 diameter_m=0.10, length_m=220.0,
                 absolute_roughness_m=roughness,
                 height_change_m=-60.0,
                 heat_transfer_coefficient_w_per_m2k=U,
                 ambient_temperature_c=T_amb,
                 n_thermal_segments=n_seg,
                 component_id="branch_a_descending"),
            Pipe(2, 4,
                 diameter_m=0.12, length_m=120.0,
                 absolute_roughness_m=roughness,
                 height_change_m=-30.0,
                 heat_transfer_coefficient_w_per_m2k=U,
                 ambient_temperature_c=T_amb,
                 n_thermal_segments=n_seg,
                 component_id="branch_b_descending"),
        ),
        node_ids=(1, 2, 3, 4),
    )
