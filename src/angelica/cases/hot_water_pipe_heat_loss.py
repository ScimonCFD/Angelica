from __future__ import annotations

from angelica.core.case import NetworkCase, PressureBoundary, ThermalBoundary
from angelica.core.components import Pipe
from angelica.properties.water_liquid import build_water_thermal_fluid


def build_hot_water_pipe_heat_loss_case() -> NetworkCase:
    """Single 500 m pipe carrying hot water (80 °C) losing heat to a 20 °C ambient.

    Network topology:
        source (node 1, 3 bar, 80 °C) ──[500 m pipe]──> sink (node 2, 1 bar)

    Fluid: water with temperature-dependent ρ(T) and μ(T).
    Pipe:  D=25 mm, U=50 W/m²K.

    The viscosity at 80 °C (~0.355 mPa·s) is roughly 3× lower than at 20 °C
    (~1.0 mPa·s), so the hydraulic solution changes as the temperature field
    evolves — requiring several outer temperature iterations to converge.

    Analytical outlet temperature (plug-flow + wall heat loss):
        T_out = T_amb + (T_in − T_amb) × exp(−U·π·D·L / (ṁ·cₚ))
    At the converged mass flow this gives roughly 39 °C.
    """
    return NetworkCase(
        name="Hot water pipe with heat loss",
        fluid_model=build_water_thermal_fluid(),
        pressure_inlets=(
            PressureBoundary(node_id=1, pressure_pa=300_000.0),
        ),
        pressure_outlets=(
            PressureBoundary(node_id=2, pressure_pa=100_000.0),
        ),
        components=(
            Pipe(
                1, 2,
                diameter_m=0.025,
                length_m=500.0,
                absolute_roughness_m=0.000045,
                heat_transfer_coefficient_w_per_m2k=50.0,
                ambient_temperature_c=20.0,
                n_thermal_segments=50,
            ),
        ),
        node_ids=(1, 2),
        thermal_inlets=(
            ThermalBoundary(node_id=1, temperature_c=80.0, bc_type="fixed_temperature"),
            ThermalBoundary(node_id=2, bc_type="zero_gradient"),
        ),
    )
