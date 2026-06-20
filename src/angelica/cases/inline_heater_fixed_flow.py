from __future__ import annotations

from angelica.core.case import FlowBoundary, NetworkCase, PressureBoundary, ThermalBoundary
from angelica.core.components import HeatSource, Pipe
from angelica.properties.thermal_fluid import ThermalFluid


def build_inline_heater_fixed_flow_case() -> NetworkCase:
    """Exact energy-balance benchmark for an inline heater with fixed mass flow.

    With adiabatic pipes, constant properties, and a prescribed mass flow, the
    heater outlet temperature should satisfy exactly:

        ΔT = Q / (ṁ cp)

    Here:
        Q = 50 kW
        ṁ = 1 kg/s
        cp = 4182 J/kg/K

    so the exact outlet temperature is 31.956 °C for a 20 °C inlet.
    """
    fluid = ThermalFluid.from_constants(
        density_kg_per_m3=998.2,
        viscosity_pa_s=1.0e-3,
        specific_heat_j_per_kg_k=4182.0,
        thermal_conductivity_w_per_m_k=0.6,
    )

    return NetworkCase(
        name="Inline heater fixed-flow benchmark",
        fluid_model=fluid,
        pressure_inlets=(),
        pressure_outlets=(
            PressureBoundary(node_id=4, pressure_pa=101_325.0),
        ),
        flow_inlets=(
            FlowBoundary(node_id=1, mass_flow_kg_per_s=1.0),
        ),
        components=(
            Pipe(
                1, 2,
                diameter_m=0.04,
                length_m=2.0,
                absolute_roughness_m=0.0,
                heat_transfer_coefficient_w_per_m2k=0.0,
                ambient_temperature_c=20.0,
                n_thermal_segments=2,
                component_id="feed_pipe",
            ),
            HeatSource(
                2, 3,
                diameter_m=0.04,
                power_w=50_000.0,
                pressure_drop_mode="rated",
                pressure_drop_pa=0.0,
                rated_mass_flow_kg_per_s=1.0,
                n_thermal_segments=10,
                length_m=1.0,
                component_id="heater",
            ),
            Pipe(
                3, 4,
                diameter_m=0.04,
                length_m=2.0,
                absolute_roughness_m=0.0,
                heat_transfer_coefficient_w_per_m2k=0.0,
                ambient_temperature_c=20.0,
                n_thermal_segments=2,
                component_id="exit_pipe",
            ),
        ),
        node_ids=(1, 2, 3, 4),
        thermal_inlets=(
            ThermalBoundary(node_id=1, temperature_c=20.0, bc_type="fixed_temperature"),
            ThermalBoundary(node_id=4, bc_type="zero_gradient"),
        ),
    )
