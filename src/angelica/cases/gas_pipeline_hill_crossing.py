from __future__ import annotations

from angelica.core.case import FlowBoundary, NetworkCase, PressureBoundary, ThermalBoundary
from angelica.core.components import Pipe
from angelica.properties.compressible_fluid import CompressibleFluid
from angelica.properties.eos import IdealGasEOS

_M_METHANE = 0.016043   # kg/mol
_MU_METHANE = 1.1e-5    # Pa·s  at 15 °C
_CP_METHANE = 2_220.0   # J/(kg·K)
_K_METHANE  = 0.033     # W/(m·K)


def build_gas_pipeline_hill_crossing_case() -> NetworkCase:
    """Natural gas pipeline crossing a significant hill — mixed elevation + outflow BC.

    The pipeline climbs from a high-pressure supply station (node 1) to a
    hilltop junction (node 2), then splits into two descending branches.
    Branch A (2→3) delivers to a fixed-pressure offtake; branch B (2→4)
    descends past the inlet elevation and exits via a prescribed mass-flow BC.

    All pipe lengths are measured along the pipe axis (actual inclined length),
    so length_m > |height_change_m| for every pipe.  The elevation correction
    ρ·g·Δz is included in every friction closure; because gas density varies
    with pressure, the gravity contribution is updated each outer iteration.

    Elevation convention (relative to node 1):
        Node 1  — 0 m      (supply station)
        Node 2  — +500 m   (hilltop junction)
        Node 3  — +200 m   (branch A offtake, 300 m below hilltop)
        Node 4  — −50 m    (branch B valley, 550 m below hilltop)

    Pipes (smooth steel, ε = 1×10⁻⁵ m):
        1 → 2  D=0.30 m  L=600 m  Δz=+500 m  (ascending trunk)
        2 → 3  D=0.25 m  L=380 m  Δz=−300 m  (branch A, partial descent)
        2 → 4  D=0.20 m  L=640 m  Δz=−550 m  (branch B, full descent past supply)

    Boundary conditions:
        Node 1 — pressure inlet  (8 000 000 Pa, T = 15 °C)
        Node 3 — pressure outlet (5 500 000 Pa)
        Node 4 — outflow BC      (ṁ_out = 8.0 kg/s, pressure free)
    """
    fluid = CompressibleFluid.from_constants(
        eos=IdealGasEOS(molecular_weight_kg_per_mol=_M_METHANE),
        viscosity_pa_s=_MU_METHANE,
        specific_heat_j_per_kg_k=_CP_METHANE,
        thermal_conductivity_w_per_m_k=_K_METHANE,
        reference_pressure_pa=7_000_000.0,
        reference_temperature_c=15.0,
    )
    roughness = 1.0e-5   # smooth steel, m

    return NetworkCase(
        name="Gas pipeline hill crossing — mixed elevation + outflow BC",
        fluid_model=fluid,
        pressure_inlets=(
            PressureBoundary(node_id=1, pressure_pa=8_000_000.0),
        ),
        pressure_outlets=(
            PressureBoundary(node_id=3, pressure_pa=5_500_000.0),
        ),
        flow_outlets=(
            FlowBoundary(node_id=4, mass_flow_kg_per_s=8.0),
        ),
        thermal_inlets=(
            ThermalBoundary(node_id=1, temperature_c=15.0),
        ),
        components=(
            Pipe(1, 2,
                 diameter_m=0.30, length_m=600.0,
                 absolute_roughness_m=roughness,
                 height_change_m=+500.0,
                 component_id="trunk_ascending"),
            Pipe(2, 3,
                 diameter_m=0.25, length_m=380.0,
                 absolute_roughness_m=roughness,
                 height_change_m=-300.0,
                 component_id="branch_a_partial_descent"),
            Pipe(2, 4,
                 diameter_m=0.20, length_m=640.0,
                 absolute_roughness_m=roughness,
                 height_change_m=-550.0,
                 component_id="branch_b_full_descent"),
        ),
        node_ids=(1, 2, 3, 4),
    )
