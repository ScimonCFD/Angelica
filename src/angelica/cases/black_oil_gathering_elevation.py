from __future__ import annotations

from angelica.core.case import NetworkCase, PressureBoundary, ThermalBoundary, InletFluidBC
from angelica.core.components import Pipe
from angelica.properties.black_oil import BlackOilFluid


def build_black_oil_gathering_elevation_case() -> NetworkCase:
    """Three-phase gathering network with wells at different topographic positions.

    Well A sits at high elevation (+150 m above the manifold) and flows
    downhill to the central manifold — gravity assists the flow.  Well B sits
    below the manifold (−100 m) and must push fluid uphill against gravity —
    requiring a higher wellhead pressure.  Both streams blend at the manifold
    and are delivered to a low-pressure separator.

    The elevation correction ρ_mix·g·Δz is significant at these height
    differences: for a 700 kg/m³ mixture it amounts to about 1 MPa per
    100 m, which visibly shifts the flow split between the two wells.

    Elevation convention (relative to manifold, node 3 = 0 m):
        Node 1 (Well A)    — +150 m  (high elevation wellhead)
        Node 2 (Well B)    — −100 m  (sub-surface wellhead)
        Node 3 (Manifold)  —    0 m  (central collection point)
        Node 4 (Separator) −  −30 m  (slight descent from manifold)

    Pipes (carbon-steel, ε = 4.6×10⁻⁵ m):
        1 → 3  D=0.15 m  L=200 m  Δz=−150 m  (well A, descends to manifold)
        2 → 3  D=0.12 m  L=180 m  Δz=+100 m  (well B, ascends to manifold)
        3 → 4  D=0.20 m  L= 80 m  Δz= −30 m  (manifold to separator)

    Fluid:
        Well A — 32°API crude,  GOR=25 m³/m³, WOR=0.4, gas_gravity=0.65
        Well B — 28°API crude,  GOR=15 m³/m³, WOR=0.8, gas_gravity=0.68

    Boundary conditions:
        Node 1 — pressure inlet  (6 000 000 Pa, T = 65 °C)
        Node 2 — pressure inlet  (9 000 000 Pa, T = 60 °C)
        Node 4 — pressure outlet (3 000 000 Pa)
    """
    roughness = 4.6e-5   # carbon-steel, m

    default_fluid = BlackOilFluid(
        api_gravity=32.0,
        gas_gravity=0.65,
        gor_sc_m3_per_m3=25.0,
        wor_sc_m3_per_m3=0.4,
        reference_pressure_pa=5_000_000.0,
        reference_temperature_c=65.0,
    )

    return NetworkCase(
        name="Black-oil gathering — mixed elevation, two wells",
        fluid_model=default_fluid,
        pressure_inlets=(
            PressureBoundary(node_id=1, pressure_pa=6_000_000.0),
            PressureBoundary(node_id=2, pressure_pa=9_000_000.0),
        ),
        pressure_outlets=(
            PressureBoundary(node_id=4, pressure_pa=3_000_000.0),
        ),
        thermal_inlets=(
            ThermalBoundary(node_id=1, temperature_c=65.0, bc_type="fixed_temperature"),
            ThermalBoundary(node_id=2, temperature_c=60.0, bc_type="fixed_temperature"),
        ),
        inlet_fluid_bcs=(
            InletFluidBC(node_id=1, api_gravity=32.0, gas_gravity=0.65,
                         gor_sc_m3_per_m3=25.0, wor_sc_m3_per_m3=0.4),
            InletFluidBC(node_id=2, api_gravity=28.0, gas_gravity=0.68,
                         gor_sc_m3_per_m3=15.0, wor_sc_m3_per_m3=0.8),
        ),
        components=(
            Pipe(1, 3,
                 diameter_m=0.15, length_m=200.0,
                 absolute_roughness_m=roughness,
                 height_change_m=-150.0,
                 heat_transfer_coefficient_w_per_m2k=5.0,
                 ambient_temperature_c=20.0,
                 n_thermal_segments=20,
                 component_id="well_a_flowline"),
            Pipe(2, 3,
                 diameter_m=0.12, length_m=180.0,
                 absolute_roughness_m=roughness,
                 height_change_m=+100.0,
                 heat_transfer_coefficient_w_per_m2k=5.0,
                 ambient_temperature_c=20.0,
                 n_thermal_segments=20,
                 component_id="well_b_flowline"),
            Pipe(3, 4,
                 diameter_m=0.20, length_m=80.0,
                 absolute_roughness_m=roughness,
                 height_change_m=-30.0,
                 heat_transfer_coefficient_w_per_m2k=5.0,
                 ambient_temperature_c=20.0,
                 n_thermal_segments=10,
                 component_id="manifold_to_separator"),
        ),
        node_ids=(1, 2, 3, 4),
    )
