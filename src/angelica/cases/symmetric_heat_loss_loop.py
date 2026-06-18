from __future__ import annotations

from angelica.core.case import FlowBoundary, NetworkCase, PressureBoundary, ThermalBoundary
from angelica.core.components import Pipe
from angelica.properties.thermal_fluid import ThermalFluid


def build_symmetric_heat_loss_loop_case() -> NetworkCase:
    """Symmetric looped network with exact branch and outlet temperature.

    Two identical heat-losing branches connect the split and merge junctions.
    By symmetry the branch mass flow rates are equal, so each branch temperature
    follows the exact NTU analytical solution and the mixed outlet temperature
    is identical to the branch outlet temperature.
    """
    fluid = ThermalFluid.from_constants(
        density_kg_per_m3=998.2,
        viscosity_pa_s=1.0e-3,
        specific_heat_j_per_kg_k=4182.0,
        thermal_conductivity_w_per_m_k=0.6,
    )

    return NetworkCase(
        name="Symmetric heat-loss loop benchmark",
        fluid_model=fluid,
        pressure_inlets=(),
        pressure_outlets=(PressureBoundary(node_id=6, pressure_pa=101_325.0),),
        flow_inlets=(FlowBoundary(node_id=1, mass_flow_kg_per_s=2.0),),
        components=(
            Pipe(1, 2, diameter_m=0.04, length_m=5.0, absolute_roughness_m=0.0,
                 heat_transfer_coefficient_w_per_m2k=0.0, ambient_temperature_c=20.0,
                 n_thermal_segments=2, component_id="inlet_header"),
            Pipe(2, 3, diameter_m=0.025, length_m=500.0, absolute_roughness_m=0.0,
                 heat_transfer_coefficient_w_per_m2k=50.0, ambient_temperature_c=20.0,
                 n_thermal_segments=120, component_id="upper_branch"),
            Pipe(2, 4, diameter_m=0.025, length_m=500.0, absolute_roughness_m=0.0,
                 heat_transfer_coefficient_w_per_m2k=50.0, ambient_temperature_c=20.0,
                 n_thermal_segments=120, component_id="lower_branch"),
            Pipe(3, 5, diameter_m=0.025, length_m=10.0, absolute_roughness_m=0.0,
                 heat_transfer_coefficient_w_per_m2k=0.0, ambient_temperature_c=20.0,
                 n_thermal_segments=2, component_id="upper_return"),
            Pipe(4, 5, diameter_m=0.025, length_m=10.0, absolute_roughness_m=0.0,
                 heat_transfer_coefficient_w_per_m2k=0.0, ambient_temperature_c=20.0,
                 n_thermal_segments=2, component_id="lower_return"),
            Pipe(5, 6, diameter_m=0.04, length_m=5.0, absolute_roughness_m=0.0,
                 heat_transfer_coefficient_w_per_m2k=0.0, ambient_temperature_c=20.0,
                 n_thermal_segments=2, component_id="outlet_header"),
        ),
        node_ids=(1, 2, 3, 4, 5, 6),
        thermal_inlets=(
            ThermalBoundary(node_id=1, temperature_c=80.0, bc_type="fixed_temperature"),
            ThermalBoundary(node_id=6, bc_type="zero_gradient"),
        ),
    )
