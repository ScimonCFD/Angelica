from __future__ import annotations

from angelica.core.case import NetworkCase, PressureBoundary, ThermalBoundary
from angelica.core.components import Pipe
from angelica.properties.water_liquid import build_water_thermal_fluid


def build_district_heating_branch_case() -> NetworkCase:
    """District heating network with one supply and two loads.

    Network topology:

        source (node 1, 4 bar, 85 °C)
             |
          [main pipe: D=50 mm, L=200 m, U=5 W/m²K]
             |
          junction (node 2)
          /         \\
     [branch A]   [branch B]
     D=25 mm      D=25 mm
     L=300 m      L=100 m
     U=5 W/m²K    U=5 W/m²K
          |               |
       sink A           sink B
      (node 3,         (node 4,
      1.5 bar)         1.5 bar)

    Fluid: water with temperature-dependent ρ(T) and μ(T) — the viscosity
    difference between 85 °C supply and ~80 °C delivery is significant enough
    to require several outer temperature iterations.

    The longer branch A loses more heat than branch B, so the delivery
    temperature at node 3 is lower than at node 4 despite both starting
    from the same junction temperature.
    """
    return NetworkCase(
        name="District heating branch",
        fluid_model=build_water_thermal_fluid(),
        pressure_inlets=(
            PressureBoundary(node_id=1, pressure_pa=400_000.0),
        ),
        pressure_outlets=(
            PressureBoundary(node_id=3, pressure_pa=150_000.0),
            PressureBoundary(node_id=4, pressure_pa=150_000.0),
        ),
        components=(
            Pipe(
                1, 2,
                diameter_m=0.05,
                length_m=200.0,
                absolute_roughness_m=0.000045,
                heat_transfer_coefficient_w_per_m2k=5.0,
                ambient_temperature_c=10.0,
                n_thermal_segments=20,
            ),
            Pipe(
                2, 3,
                diameter_m=0.025,
                length_m=300.0,
                absolute_roughness_m=0.000045,
                heat_transfer_coefficient_w_per_m2k=5.0,
                ambient_temperature_c=10.0,
                n_thermal_segments=20,
            ),
            Pipe(
                2, 4,
                diameter_m=0.025,
                length_m=100.0,
                absolute_roughness_m=0.000045,
                heat_transfer_coefficient_w_per_m2k=5.0,
                ambient_temperature_c=10.0,
                n_thermal_segments=20,
            ),
        ),
        node_ids=(1, 2, 3, 4),
        thermal_inlets=(
            ThermalBoundary(node_id=1, temperature_c=85.0, bc_type="fixed_temperature"),
            ThermalBoundary(node_id=3, bc_type="zero_gradient"),
            ThermalBoundary(node_id=4, bc_type="zero_gradient"),
        ),
    )
