from __future__ import annotations

from angelica.core.case import NetworkCase, PressureBoundary
from angelica.core.components import Fitting, Pipe
from angelica.properties.single_component import SingleComponentFluid


def build_steady_water_network_case() -> NetworkCase:
    return NetworkCase(
        name="Recovered steady water network",
        fluid_model=SingleComponentFluid(
            density_kg_per_m3=998.25,
            viscosity_pa_s=0.001,
        ),
        pressure_inlets=(
            PressureBoundary(node_id=1, pressure_pa=251.3 * 1000.0),
            PressureBoundary(node_id=5, pressure_pa=201.3 * 1000.0),
        ),
        pressure_outlets=(
            PressureBoundary(node_id=6, pressure_pa=101.3 * 1000.0),
        ),
        components=(
            Pipe(1, 7, 0.05, length_m=150.0, absolute_roughness_m=0.000045, height_change_m=0.0),
            Fitting(7, 8, 0.05, loss_coefficient=1.5),
            Pipe(8, 2, 0.05, length_m=50.0, absolute_roughness_m=0.000045, height_change_m=0.0),
            Pipe(2, 9, 0.025, length_m=10.0, absolute_roughness_m=0.000045, height_change_m=0.0),
            Fitting(9, 10, 0.025, loss_coefficient=1.5),
            Pipe(10, 3, 0.025, length_m=20.0, absolute_roughness_m=0.000045, height_change_m=0.0),
            Pipe(2, 4, 0.025, length_m=40.0, absolute_roughness_m=0.000045, height_change_m=0.0),
            Pipe(4, 3, 0.025, length_m=15.0, absolute_roughness_m=0.000045, height_change_m=0.0),
            Pipe(5, 4, 0.05, length_m=300.0, absolute_roughness_m=0.000045, height_change_m=0.0),
            Pipe(3, 6, 0.05, length_m=60.0, absolute_roughness_m=0.000045, height_change_m=0.0),
        ),
        node_ids=(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
    )
