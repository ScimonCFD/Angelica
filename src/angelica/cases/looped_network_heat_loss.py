from __future__ import annotations

from angelica.core.case import NetworkCase, PressureBoundary, ThermalBoundary
from angelica.core.components import Pipe
from angelica.properties.water_liquid import build_water_thermal_fluid


def build_looped_network_heat_loss_case() -> NetworkCase:
    """Looped pipe network with heat loss — illustrates outer-loop coupling.

    Network topology:

        source (node 1, 5 bar, 95 °C)
             |
          [feed: D=50 mm, L=100 m]
             |
          junction A (node 2)
          /                   \\
     [upper branch]       [direct bypass]
     2→3 D=40mm 400m       2→4 D=35mm 400m
     3→4 D=40mm 300m
          \\                  /
          junction C (node 4)
             |
          [exit: D=50 mm, L=100 m]
             |
          sink (node 5, 1 bar)

    All pipes: U = 30 W/m²K, T_amb = 10 °C.

    The upper branch (700 m total) loses more heat than the direct bypass
    (400 m), creating different average temperatures and viscosities on each
    path.  This viscosity contrast couples the hydraulic and thermal solutions,
    so the outer temperature loop requires multiple iterations even with a
    near-balanced network.
    """
    return NetworkCase(
        name="Looped network with heat loss",
        fluid_model=build_water_thermal_fluid(),
        pressure_inlets=(
            PressureBoundary(node_id=1, pressure_pa=500_000.0),
        ),
        pressure_outlets=(
            PressureBoundary(node_id=5, pressure_pa=100_000.0),
        ),
        components=(
            Pipe(1, 2, diameter_m=0.05, length_m=100.0,
                 absolute_roughness_m=0.000045,
                 heat_transfer_coefficient_w_per_m2k=30.0,
                 ambient_temperature_c=10.0, n_thermal_segments=1),
            Pipe(2, 3, diameter_m=0.04, length_m=400.0,
                 absolute_roughness_m=0.000045,
                 heat_transfer_coefficient_w_per_m2k=30.0,
                 ambient_temperature_c=10.0, n_thermal_segments=1),
            Pipe(3, 4, diameter_m=0.04, length_m=300.0,
                 absolute_roughness_m=0.000045,
                 heat_transfer_coefficient_w_per_m2k=30.0,
                 ambient_temperature_c=10.0, n_thermal_segments=1),
            Pipe(2, 4, diameter_m=0.035, length_m=400.0,
                 absolute_roughness_m=0.000045,
                 heat_transfer_coefficient_w_per_m2k=30.0,
                 ambient_temperature_c=10.0, n_thermal_segments=1),
            Pipe(4, 5, diameter_m=0.05, length_m=100.0,
                 absolute_roughness_m=0.000045,
                 heat_transfer_coefficient_w_per_m2k=30.0,
                 ambient_temperature_c=10.0, n_thermal_segments=1),
        ),
        node_ids=(1, 2, 3, 4, 5),
        thermal_inlets=(
            ThermalBoundary(node_id=1, temperature_c=95.0, bc_type="fixed_temperature"),
            ThermalBoundary(node_id=5, bc_type="zero_gradient"),
        ),
    )
