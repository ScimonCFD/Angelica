"""Tutorial 07 — Composite U loop: insulated vs bare arm heat loss.

Demonstrates the multi-layer thermal resistance model (composite U) on a
looped crude-oil pipeline.  The loop splits into two parallel arms:

  - Upper arm (2→3→4): mineral-wool–insulated steel pipe
      U_overall ≈ 0.74 W/m²K
  - Lower arm (2→4):   bare uninsulated steel pipe
      U_overall ≈ 14.5 W/m²K

Both arms start at the same temperature (80 °C) and arrive at the merge node
(Node 4) at significantly different exit temperatures, illustrating how
insulation quality governs heat loss in looped networks.

Topology
--------
  1 ─── A ─── 2 ─── B ─── 3 ─── C ─── 4 ─── E ─── 5
                └────────── D ───────────┘
"""
from __future__ import annotations

import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[4]  # …/angelica/angelica/src
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.core.case import NetworkCase, PressureBoundary, ThermalBoundary
from angelica.core.components import Pipe
from angelica.core.settings import SolverSettings
from angelica.properties.thermal_fluid import ThermalFluid
from angelica.solvers import NonIsothermalSolverSettings, SteadyNonIsothermalIncompressibleSolver

# ── Fluid (crude oil at ~80 °C) ───────────────────────────────────────────────
FLUID = ThermalFluid.from_constants(
    density_kg_per_m3=855.0,
    viscosity_pa_s=6e-3,
    specific_heat_j_per_kg_k=2000.0,
    thermal_conductivity_w_per_m_k=0.13,
)

T_INLET = 80.0   # °C — hot crude arriving from upstream
T_AMB   = 10.0   # °C — cold ambient (buried / exposed)

# ── Composite layer parameters ────────────────────────────────────────────────
# Insulated: carbon-steel pipe + 50 mm mineral wool + cladding
# R = 1/500 + 0.008/50 + 0.05/0.04 + 1/10 ≈ 1.352 m²K/W  → U ≈ 0.74 W/m²K
_INSULATED = dict(
    inner_film_coefficient_w_per_m2k=500.0,
    wall_thickness_m=0.008,
    wall_thermal_conductivity_w_per_mk=50.0,
    insulation_thickness_m=0.050,
    insulation_thermal_conductivity_w_per_mk=0.04,
    outer_film_coefficient_w_per_m2k=10.0,
    ambient_temperature_c=T_AMB,
)

# Bare: carbon-steel pipe with outer wind-exposed film only
# R = 1/500 + 0.008/50 + 1/15 ≈ 0.069 m²K/W  → U ≈ 14.5 W/m²K
_BARE = dict(
    inner_film_coefficient_w_per_m2k=500.0,
    wall_thickness_m=0.008,
    wall_thermal_conductivity_w_per_mk=50.0,
    insulation_thickness_m=0.0,
    insulation_thermal_conductivity_w_per_mk=0.0,
    outer_film_coefficient_w_per_m2k=15.0,
    ambient_temperature_c=T_AMB,
)


def build_case() -> NetworkCase:
    pipes = [
        Pipe(start_node=1, end_node=2, diameter_m=0.10, length_m=200.0,
             absolute_roughness_m=4.6e-5,
             heat_transfer_coefficient_w_per_m2k=5.0, ambient_temperature_c=T_AMB,
             n_thermal_segments=5, component_id="A"),
        # Upper arm — insulated (composite U ≈ 0.74 W/m²K)
        Pipe(start_node=2, end_node=3, diameter_m=0.08, length_m=800.0,
             absolute_roughness_m=4.6e-5, n_thermal_segments=10,
             component_id="B", **_INSULATED),
        Pipe(start_node=3, end_node=4, diameter_m=0.08, length_m=800.0,
             absolute_roughness_m=4.6e-5, n_thermal_segments=10,
             component_id="C", **_INSULATED),
        # Lower arm — bare steel (composite U ≈ 14.5 W/m²K)
        Pipe(start_node=2, end_node=4, diameter_m=0.10, length_m=1200.0,
             absolute_roughness_m=4.6e-5, n_thermal_segments=15,
             component_id="D", **_BARE),
        Pipe(start_node=4, end_node=5, diameter_m=0.10, length_m=200.0,
             absolute_roughness_m=4.6e-5,
             heat_transfer_coefficient_w_per_m2k=5.0, ambient_temperature_c=T_AMB,
             n_thermal_segments=5, component_id="E"),
    ]

    return NetworkCase(
        name="Composite-U loop — insulated vs bare arm",
        fluid_model=FLUID,
        pressure_inlets=(PressureBoundary(node_id=1, pressure_pa=12e5),),
        pressure_outlets=(PressureBoundary(node_id=5, pressure_pa=2e5),),
        components=tuple(pipes),
        node_ids=(1, 2, 3, 4, 5),
        thermal_inlets=(
            ThermalBoundary(node_id=1, bc_type="fixed_temperature", temperature_c=T_INLET),
            ThermalBoundary(node_id=5, bc_type="zero_gradient"),
        ),
    )


def main() -> None:
    case = build_case()
    solver = SteadyNonIsothermalIncompressibleSolver(
        hydraulic_settings=SolverSettings(
            turbulent_iterations=150,
            pressure_correction_abs_tolerance_pa=1e-2,
        ),
        non_isothermal_settings=NonIsothermalSolverSettings(
            max_temperature_iterations=60,
            temperature_tolerance_k=0.01,
        ),
    )
    result = solver.solve(case)

    print(f"Case: {case.name}")
    print(f"Converged: {result.converged}")
    print()

    node_labels = {
        1: "Source",
        2: "Split junction",
        3: "Upper-arm midpoint",
        4: "Merge junction",
        5: "Sink",
    }
    print(f"{'Node':<26}  {'P (bar)':>8}  {'T (°C)':>8}")
    for nid in sorted(result.node_pressures_pa):
        p = result.node_pressures_pa[nid] / 1e5
        T = result.node_temperatures_c[nid]
        print(f"  {node_labels[nid]:<24}  {p:>8.3f}  {T:>8.2f}")
    print()

    def composite_u(inner_h, t_wall, k_wall, t_ins, k_ins, outer_h):
        R = 0.0
        if inner_h > 0: R += 1.0 / inner_h
        if t_wall > 0 and k_wall > 0: R += t_wall / k_wall
        if t_ins > 0 and k_ins > 0: R += t_ins / k_ins
        if outer_h > 0: R += 1.0 / outer_h
        return 1.0 / R if R > 0 else 0.0

    U_ins  = composite_u(500, 0.008, 50.0, 0.05, 0.04, 10.0)
    U_bare = composite_u(500, 0.008, 50.0, 0.00, 0.00, 15.0)
    print(f"Effective U — insulated arm : {U_ins:.2f} W/m²K")
    print(f"Effective U — bare arm      : {U_bare:.2f} W/m²K")
    print(f"Ratio (bare/insulated)      : {U_bare/U_ins:.1f}×")
    print()

    pipe_labels = {
        "A": "Feed A  1→2 (D100, 200m, U=5.0)",
        "B": "Upper B 2→3 (D80, 800m, insul)",
        "C": "Upper C 3→4 (D80, 800m, insul)",
        "D": "Bare D  2→4 (D100, 1200m, bare)",
        "E": "Exit E  4→5 (D100, 200m, U=5.0)",
    }
    header = f"{'Pipe':<34}  {'kg/s':>8}  {'m³/h':>8}  {'T_in °C':>8}  {'T_out °C':>8}"
    print(header)
    for cf in result.component_flows:
        pid = cf.label.split(":")[-1]
        label = pipe_labels.get(pid, cf.label)
        T_in  = f"{cf.temperature_in_c:.1f}"  if cf.temperature_in_c  is not None else "—"
        T_out = f"{cf.temperature_out_c:.1f}" if cf.temperature_out_c is not None else "—"
        print(f"  {label:<32}  {cf.mass_flow_kg_per_s:>8.4f}  {cf.volumetric_flow_m3_per_h:>8.3f}"
              f"  {T_in:>8}  {T_out:>8}")


if __name__ == "__main__":
    main()
