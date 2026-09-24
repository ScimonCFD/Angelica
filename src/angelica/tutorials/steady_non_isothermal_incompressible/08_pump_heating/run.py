"""Tutorial 08 — Pump shaft-work heating in a crude oil pipeline.

A crude oil gathering system collects oil at 70 °C and feeds a pump station
that boosts the oil into a high-pressure export line.  This tutorial shows:

  1. The pump raises fluid temperature by ΔT = ΔP / (ρ · Cp).
  2. The pump heating partially offsets the heat lost in the upstream pipe.
  3. The analytical shaft-work formula matches the solver result exactly.

Network layout (series):

    Source ──── Pipe 1 ──── Pump ──── Pipe 2 ──── Sink
    (2 bar,         500 m,                 1000 m,   12 bar)
    70 °C        U=3 W/m²K            U=3 W/m²K
"""
from __future__ import annotations

import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[4]  # …/angelica/angelica/src
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.core.case import NetworkCase, PressureBoundary, ThermalBoundary
from angelica.core.components import Pipe, Pump
from angelica.properties.thermal_fluid import ThermalFluid
from angelica.solvers.steady_non_isothermal_incompressible import SteadyNonIsothermalIncompressibleSolver

# ── Fluid: medium crude oil at ~70 °C ────────────────────────────────────────
_RHO = 870.0    # kg/m³
_CP  = 2000.0   # J/(kg·K)

_FLUID = ThermalFluid.from_constants(
    density_kg_per_m3=_RHO,
    viscosity_pa_s=6e-3,
    specific_heat_j_per_kg_k=_CP,
    thermal_conductivity_w_per_m_k=0.13,
)

_D        = 0.15    # m  — pipe diameter
_U        = 3.0     # W/(m²·K) — light insulation (buried pipe)
_T_AMB    = 15.0    # °C — soil temperature
_T_INLET  = 70.0    # °C — gathering temperature
_P_SOURCE = 2e5     # Pa — inlet pressure
_P_SINK   = 12e5    # Pa — export pipeline pressure


def build_case() -> NetworkCase:
    return NetworkCase(
        name="Pump shaft-work heating — crude oil",
        fluid_model=_FLUID,
        node_ids=(1, 2, 3, 4),
        pressure_inlets=(PressureBoundary(node_id=1, pressure_pa=_P_SOURCE),),
        pressure_outlets=(PressureBoundary(node_id=4, pressure_pa=_P_SINK),),
        thermal_inlets=(
            ThermalBoundary(node_id=1, bc_type="fixed_temperature", temperature_c=_T_INLET),
            ThermalBoundary(node_id=4, bc_type="zero_gradient"),
        ),
        components=(
            Pipe(
                start_node=1, end_node=2,
                diameter_m=_D, length_m=500.0,
                absolute_roughness_m=4.6e-5,
                heat_transfer_coefficient_w_per_m2k=_U,
                ambient_temperature_c=_T_AMB,
                n_thermal_segments=10,
                component_id="upstream_pipe",
            ),
            Pump(
                start_node=2, end_node=3,
                diameter_m=_D,
                curve_points_q_head=(
                    (0.000, 130.0),
                    (0.005, 129.0),
                    (0.010, 127.0),
                    (0.020, 122.0),
                    (0.030, 115.0),
                    (0.040,  95.0),
                    (0.050,   0.0),
                ),
                component_id="pump",
            ),
            Pipe(
                start_node=3, end_node=4,
                diameter_m=_D, length_m=1000.0,
                absolute_roughness_m=4.6e-5,
                heat_transfer_coefficient_w_per_m2k=_U,
                ambient_temperature_c=_T_AMB,
                n_thermal_segments=20,
                component_id="downstream_pipe",
            ),
        ),
    )


def main() -> None:
    result = SteadyNonIsothermalIncompressibleSolver().solve(build_case())

    print(f"Case: {result.converged and 'Converged' or 'NOT CONVERGED'}")
    print()

    T1 = result.node_temperatures_c[1]
    T2 = result.node_temperatures_c[2]
    T3 = result.node_temperatures_c[3]
    T4 = result.node_temperatures_c[4]
    P2 = result.node_pressures_pa[2]
    P3 = result.node_pressures_pa[3]

    mdot = next(f.mass_flow_kg_per_s for f in result.component_flows if f.label and "pump" in f.label.lower())

    print("Node temperatures:")
    print(f"  Node 1 — source inlet:       {T1:.3f} °C  ({result.node_pressures_pa[1]/1e5:.2f} bar)")
    print(f"  Node 2 — upstream pump:      {T2:.3f} °C  ({P2/1e5:.2f} bar)")
    print(f"  Node 3 — downstream pump:    {T3:.3f} °C  ({P3/1e5:.2f} bar)")
    print(f"  Node 4 — export sink:        {T4:.3f} °C  ({result.node_pressures_pa[4]/1e5:.2f} bar)")
    print()

    dT_pipe1 = T2 - T1
    dT_pump  = T3 - T2
    dT_pipe2 = T4 - T3
    print("Temperature changes:")
    print(f"  Pipe 1 (500 m heat loss):   {dT_pipe1:+.3f} °C")
    print(f"  Pump (shaft work):          {dT_pump:+.3f} °C")
    print(f"  Pipe 2 (1000 m heat loss):  {dT_pipe2:+.3f} °C")
    print()

    dP_pump = P3 - P2
    dT_analytical = dP_pump / (_RHO * _CP)
    print("Shaft-work verification:")
    print(f"  ΔP pump:                    {dP_pump/1e5:.3f} bar")
    print(f"  ṁ:                          {mdot:.4f} kg/s")
    print(f"  ΔT analytical = ΔP/(ρ·Cp): {dT_analytical:+.4f} °C")
    print(f"  ΔT solver:                  {dT_pump:+.4f} °C")
    print(f"  Difference:                 {abs(dT_pump - dT_analytical):.2e} °C")


if __name__ == "__main__":
    main()
