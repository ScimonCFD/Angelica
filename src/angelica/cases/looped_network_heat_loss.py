from __future__ import annotations

from angelica.core.case import NetworkCase, PressureBoundary, ThermalBoundary
from angelica.core.components import Pipe
from angelica.properties.water_liquid import build_water_thermal_fluid


def build_looped_network_heat_loss_case() -> NetworkCase:
    """Looped pipe network with significant heat loss — designed to require
    multiple outer temperature iterations.

    Network topology:

        source (node 1, 6 bar, 95 °C)
             |
          [feed pipe: D=50 mm, L=50 m]
             |
          junction A (node 2)
          /                    \\
     [long path]           [short bypass]
     2→3 D=20mm 600m        2→4 D=25mm 100m
     3→4 D=20mm 400m
          \\                  /
          junction C (node 4)    ← hot bypass + cold long path mix here
             |
          [exit pipe: D=50 mm, L=50 m]
             |
          sink (node 5, 1 bar)

    The long path (1000 m total, D=20 mm) loses almost all its heat to the
    5 °C ambient — arriving at junction C near-cold (~21 °C).  The short
    bypass (100 m, D=25 mm) loses very little heat and arrives hot (~90 °C).
    The strong viscosity difference between the two paths (μ is ~2× higher
    in the cold long path) creates a tight coupling between the hydraulic
    and thermal fields that forces the outer temperature loop to iterate.

    With temperature_relaxation=0.5, the solver takes ~14 outer iterations
    to converge, producing a clear, smooth convergence curve.
    """
    return NetworkCase(
        name="Looped network with heat loss",
        fluid_model=build_water_thermal_fluid(),
        pressure_inlets=(
            PressureBoundary(node_id=1, pressure_pa=600_000.0),
        ),
        pressure_outlets=(
            PressureBoundary(node_id=5, pressure_pa=100_000.0),
        ),
        components=(
            Pipe(1, 2, diameter_m=0.05,  length_m=50.0,
                 absolute_roughness_m=0.000045,
                 heat_transfer_coefficient_w_per_m2k=50.0,
                 ambient_temperature_c=5.0, n_thermal_segments=10),
            Pipe(2, 3, diameter_m=0.02,  length_m=600.0,
                 absolute_roughness_m=0.000045,
                 heat_transfer_coefficient_w_per_m2k=50.0,
                 ambient_temperature_c=5.0, n_thermal_segments=30),
            Pipe(3, 4, diameter_m=0.02,  length_m=400.0,
                 absolute_roughness_m=0.000045,
                 heat_transfer_coefficient_w_per_m2k=50.0,
                 ambient_temperature_c=5.0, n_thermal_segments=20),
            Pipe(2, 4, diameter_m=0.025, length_m=100.0,
                 absolute_roughness_m=0.000045,
                 heat_transfer_coefficient_w_per_m2k=50.0,
                 ambient_temperature_c=5.0, n_thermal_segments=10),
            Pipe(4, 5, diameter_m=0.05,  length_m=50.0,
                 absolute_roughness_m=0.000045,
                 heat_transfer_coefficient_w_per_m2k=50.0,
                 ambient_temperature_c=5.0, n_thermal_segments=10),
        ),
        node_ids=(1, 2, 3, 4, 5),
        thermal_inlets=(
            ThermalBoundary(node_id=1, temperature_c=95.0, bc_type="fixed_temperature"),
            ThermalBoundary(node_id=5, bc_type="zero_gradient"),
        ),
    )
