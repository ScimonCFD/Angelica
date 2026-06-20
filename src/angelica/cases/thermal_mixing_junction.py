from __future__ import annotations

from angelica.core.case import FlowBoundary, NetworkCase, PressureBoundary, ThermalBoundary
from angelica.core.components import Pipe
from angelica.properties.thermal_fluid import ThermalFluid


def build_thermal_mixing_junction_case() -> NetworkCase:
    """Exact adiabatic mixing benchmark for the non-isothermal solver.

    Two inlet streams with prescribed mass flow rates and temperatures merge at
    a junction and leave through a single adiabatic outlet pipe.

    Topology:
        node 1 (ṁ = 1 kg/s, 80 °C) ──┐
                                      ├── node 3 ── node 4 (pressure outlet)
        node 2 (ṁ = 2 kg/s, 20 °C) ──┘

    With adiabatic pipes and constant cp, the exact mixed temperature is:

        T_mix = (ṁ1 T1 + ṁ2 T2) / (ṁ1 + ṁ2) = 40 °C
    """
    fluid = ThermalFluid.from_constants(
        density_kg_per_m3=998.2,
        viscosity_pa_s=1.0e-3,
        specific_heat_j_per_kg_k=4182.0,
        thermal_conductivity_w_per_m_k=0.6,
    )

    return NetworkCase(
        name="Thermal mixing junction benchmark",
        fluid_model=fluid,
        pressure_inlets=(),
        pressure_outlets=(
            PressureBoundary(node_id=4, pressure_pa=101_325.0),
        ),
        flow_inlets=(
            FlowBoundary(node_id=1, mass_flow_kg_per_s=1.0),
            FlowBoundary(node_id=2, mass_flow_kg_per_s=2.0),
        ),
        components=(
            Pipe(
                1, 3,
                diameter_m=0.03,
                length_m=5.0,
                absolute_roughness_m=0.0,
                heat_transfer_coefficient_w_per_m2k=0.0,
                ambient_temperature_c=20.0,
                n_thermal_segments=4,
                component_id="hot_branch",
            ),
            Pipe(
                2, 3,
                diameter_m=0.03,
                length_m=5.0,
                absolute_roughness_m=0.0,
                heat_transfer_coefficient_w_per_m2k=0.0,
                ambient_temperature_c=20.0,
                n_thermal_segments=4,
                component_id="cold_branch",
            ),
            Pipe(
                3, 4,
                diameter_m=0.04,
                length_m=5.0,
                absolute_roughness_m=0.0,
                heat_transfer_coefficient_w_per_m2k=0.0,
                ambient_temperature_c=20.0,
                n_thermal_segments=4,
                component_id="mixed_outlet",
            ),
        ),
        node_ids=(1, 2, 3, 4),
        thermal_inlets=(
            ThermalBoundary(node_id=1, temperature_c=80.0, bc_type="fixed_temperature"),
            ThermalBoundary(node_id=2, temperature_c=20.0, bc_type="fixed_temperature"),
            ThermalBoundary(node_id=4, bc_type="zero_gradient"),
        ),
    )
