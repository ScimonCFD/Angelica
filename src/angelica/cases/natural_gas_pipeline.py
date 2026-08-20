from __future__ import annotations

from angelica.core.case import NetworkCase, PressureBoundary, ThermalBoundary
from angelica.core.components import Pipe
from angelica.properties.compressible_fluid import CompressibleFluid
from angelica.properties.eos import IdealGasEOS

_M_METHANE = 0.016043       # kg/mol — molecular weight of CH₄
_MU_METHANE = 1.1e-5        # Pa·s  — dynamic viscosity at 15 °C
_CP_METHANE = 2_220.0       # J/(kg·K) — specific heat at 15 °C
_K_METHANE = 0.033          # W/(m·K) — thermal conductivity at 15 °C
_ROUGHNESS = 4.6e-5         # m — commercial steel pipe


def build_natural_gas_pipeline_case() -> NetworkCase:
    """Branched natural gas (methane) gathering network — 800 kPa inlet, 500 kPa outlets.

    Demonstrates compressible flow: inlet density ≈ 5.26 kg/m³,
    outlet density ≈ 3.29 kg/m³  (ideal gas, T = 15 °C).

    Node IDs:
        1 — Pressure inlet       (800 000 Pa)
        2 — Interior junction
        3 — Pressure outlet A    (500 000 Pa)
        4 — Pressure outlet B    (500 000 Pa)

    Pipes (commercial steel, ε = 4.6 × 10⁻⁵ m):
        Pipe 1:  1 → 2,  D = 0.15 m, L = 500 m
        Pipe 2:  2 → 3,  D = 0.10 m, L = 300 m
        Pipe 3:  2 → 4,  D = 0.10 m, L = 200 m
    """
    fluid = CompressibleFluid.from_constants(
        eos=IdealGasEOS(molecular_weight_kg_per_mol=_M_METHANE),
        viscosity_pa_s=_MU_METHANE,
        specific_heat_j_per_kg_k=_CP_METHANE,
        thermal_conductivity_w_per_m_k=_K_METHANE,
        reference_pressure_pa=500_000.0,
        reference_temperature_c=15.0,
    )

    return NetworkCase(
        name="Natural Gas Pipeline (methane, 800→500 kPa)",
        fluid_model=fluid,
        pressure_inlets=(
            PressureBoundary(node_id=1, pressure_pa=800_000.0),
        ),
        thermal_inlets=(
            ThermalBoundary(node_id=1, temperature_c=15.0),
        ),
        pressure_outlets=(
            PressureBoundary(node_id=3, pressure_pa=500_000.0),
            PressureBoundary(node_id=4, pressure_pa=500_000.0),
        ),
        components=(
            Pipe(1, 2, diameter_m=0.15, length_m=500.0,
                 absolute_roughness_m=_ROUGHNESS, component_id="pipe_1"),
            Pipe(2, 3, diameter_m=0.10, length_m=300.0,
                 absolute_roughness_m=_ROUGHNESS, component_id="pipe_2"),
            Pipe(2, 4, diameter_m=0.10, length_m=200.0,
                 absolute_roughness_m=_ROUGHNESS, component_id="pipe_3"),
        ),
        node_ids=(1, 2, 3, 4),
    )
