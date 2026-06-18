from __future__ import annotations

from angelica.core.case import NetworkCase, PressureBoundary, ThermalBoundary
from angelica.core.components import HeatSource, Pipe
from angelica.properties.water_liquid import build_water_thermal_fluid


def build_inline_heater_case() -> NetworkCase:
    """Cold water heated by an inline electric heater.

    Network topology:
        source (node 1, 3 bar, 20 °C)
            ──[feed pipe, 10 m]──
        node 2 (junction)
            ──[heater, Q=50 kW, ΔP=0]──
        node 3 (junction)
            ──[exit pipe, 10 m]──
        sink (node 4, 1 bar)

    Fluid: water with temperature-dependent ρ(T) and μ(T).
    Heater: power_w = 50 000 W, hydraulically transparent (pressure_drop_pa = 0).

    Expected temperature rise: ΔT = Q / (ṁ·cₚ).
    At ṁ ≈ 0.8 kg/s and cₚ = 4182 J/kg·K → ΔT ≈ 15 K → T_out ≈ 35 °C.
    """
    return NetworkCase(
        name="Inline heater",
        fluid_model=build_water_thermal_fluid(),
        pressure_inlets=(PressureBoundary(node_id=1, pressure_pa=300_000.0),),
        pressure_outlets=(PressureBoundary(node_id=4, pressure_pa=100_000.0),),
        components=(
            Pipe(
                1, 2,
                diameter_m=0.05,
                length_m=10.0,
                absolute_roughness_m=0.000045,
                heat_transfer_coefficient_w_per_m2k=0.0,
                ambient_temperature_c=20.0,
                n_thermal_segments=2,
            ),
            HeatSource(
                2, 3,
                diameter_m=0.05,
                power_w=50_000.0,
                pressure_drop_mode="rated",
                pressure_drop_pa=0.0,
                rated_mass_flow_kg_per_s=1.0,
                n_thermal_segments=10,
            ),
            Pipe(
                3, 4,
                diameter_m=0.05,
                length_m=10.0,
                absolute_roughness_m=0.000045,
                heat_transfer_coefficient_w_per_m2k=0.0,
                ambient_temperature_c=20.0,
                n_thermal_segments=2,
            ),
        ),
        node_ids=(1, 2, 3, 4),
        thermal_inlets=(
            ThermalBoundary(node_id=1, temperature_c=20.0, bc_type="fixed_temperature"),
            ThermalBoundary(node_id=4, bc_type="zero_gradient"),
        ),
    )
