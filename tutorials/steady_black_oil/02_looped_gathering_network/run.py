"""
Tutorial: Black-Oil Looped Gathering Network
=============================================
Two parallel flow paths connect a high-pressure inlet to a low-pressure
separator, forming a loop.  The solver distributes flow between the paths
according to their hydraulic resistances, and the PVT state differs in
each branch because the intermediate pressures are different.

Geometry
--------
                PipeA (5 km, D=0.20 m)
  [Node 1] ─────────────────────────────> [Node 2] ─┐
      │                                              │
      │  PipeC (8 km, D=0.15 m)                  PipeB
      └────────────────────────> [Node 3]          (5 km,
                                    │             D=0.18 m)
                                 PipeD              │
                              (3 km, D=0.20 m)      │
                                    └───────────────> [Node 4]
                                                   Separator

  Node 1 : P = 8 MPa, T = 60 °C  (inlet)
  Node 4 : P = 2 MPa              (outlet)
  Nodes 2, 3 : free junctions

Upper path  1→2→4 : 10 km total, wider diameters
Lower path  1→3→4 : 11 km total, narrower middle section

Fluid
-----
  32 °API crude oil, γ_g = 0.65, GOR = 25 m³/m³, WOR = 0.5 m³/m³
  Bubble point at 60 °C ≈ 5.6 MPa.  The intermediate node pressures
  determine whether each branch is undersaturated or two-phase.
"""
import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[4] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.core.case import NetworkCase, PressureBoundary, ThermalBoundary
from angelica.core.components import Pipe
from angelica.properties.black_oil import BlackOilFluid, bubble_point_pa
from angelica.solvers import SteadyBlackOilSolver

# ── Fluid ─────────────────────────────────────────────────────────────────────
API       = 32.0
GAS_GR    = 0.65
GOR_SC    = 25.0
WOR_SC    = 0.5
T_INLET_C = 60.0

fluid = BlackOilFluid(
    api_gravity             = API,
    gas_gravity             = GAS_GR,
    gor_sc_m3_per_m3        = GOR_SC,
    wor_sc_m3_per_m3        = WOR_SC,
    reference_pressure_pa   = 5e6,
    reference_temperature_c = T_INLET_C,
)

Pb = bubble_point_pa(GOR_SC, GAS_GR, API, T_INLET_C)
print(f"Bubble point at {T_INLET_C:.0f} °C : {Pb/1e6:.2f} MPa")
print()

# ── Network ───────────────────────────────────────────────────────────────────
pipes = [
    Pipe(component_id="pipe_A_1_2", start_node=1, end_node=2,
         diameter_m=0.20, length_m=5_000.0, absolute_roughness_m=46e-6,
         heat_transfer_coefficient_w_per_m2k=5.0, ambient_temperature_c=15.0),

    Pipe(component_id="pipe_B_2_4", start_node=2, end_node=4,
         diameter_m=0.18, length_m=5_000.0, absolute_roughness_m=46e-6,
         heat_transfer_coefficient_w_per_m2k=5.0, ambient_temperature_c=15.0),

    Pipe(component_id="pipe_C_1_3", start_node=1, end_node=3,
         diameter_m=0.15, length_m=8_000.0, absolute_roughness_m=46e-6,
         heat_transfer_coefficient_w_per_m2k=5.0, ambient_temperature_c=15.0),

    Pipe(component_id="pipe_D_3_4", start_node=3, end_node=4,
         diameter_m=0.20, length_m=3_000.0, absolute_roughness_m=46e-6,
         heat_transfer_coefficient_w_per_m2k=5.0, ambient_temperature_c=15.0),
]

case = NetworkCase(
    name             = "Black-Oil Looped Gathering Network",
    fluid_model      = fluid,
    pressure_inlets  = (PressureBoundary(node_id=1, pressure_pa=8e6),),
    pressure_outlets = (PressureBoundary(node_id=4, pressure_pa=2e6),),
    components       = tuple(pipes),
    thermal_inlets   = (
        ThermalBoundary(node_id=1, temperature_c=T_INLET_C, bc_type="fixed_temperature"),
    ),
)

# ── Solve ─────────────────────────────────────────────────────────────────────
solver = SteadyBlackOilSolver()
result = solver.solve(case)

print("=" * 60)
print(f"Case: {case.name}")
print(f"Converged:      {result.converged}")
print(f"PVT iterations: {len(result.density_history)}")
print()

# ── Node pressures and temperatures ──────────────────────────────────────────
print(f"{'Node':>6}  {'P (MPa)':>9}  {'T (°C)':>8}  {'Regime':>14}")
for node_id in sorted(result.node_pressures_pa):
    P = result.node_pressures_pa[node_id]
    T = result.node_temperatures_c[node_id]
    regime = "UNDERSATURATED" if P >= Pb else "two-phase"
    print(f"{node_id:>6}  {P/1e6:>9.3f}  {T:>8.2f}  {regime:>14}")
print()

# ── Per-pipe flows ────────────────────────────────────────────────────────────
print(f"{'Pipe':30s}  {'ṁ (kg/s)':>10}  {'Q_mix (m³/h)':>13}")
upper_mdot = lower_mdot = 0.0
for cf in result.component_flows:
    print(f"{cf.label:30s}  {cf.mass_flow_kg_per_s:>10.3f}  {cf.volumetric_flow_m3_per_h:>13.2f}")
    if "A" in cf.label or "B" in cf.label:
        upper_mdot = max(upper_mdot, cf.mass_flow_kg_per_s)
    else:
        lower_mdot = max(lower_mdot, cf.mass_flow_kg_per_s)

total_mdot = upper_mdot + lower_mdot
print()
print(f"Upper path (A+B): {upper_mdot:.3f} kg/s  ({100*upper_mdot/total_mdot:.1f} %)")
print(f"Lower path (C+D): {lower_mdot:.3f} kg/s  ({100*lower_mdot/total_mdot:.1f} %)")
print(f"Total:            {total_mdot:.3f} kg/s")

# ── Surface rates (total) ─────────────────────────────────────────────────────
from angelica.properties.dead_oil import dead_oil_density_kg_per_m3

rho_oil_sc = dead_oil_density_kg_per_m3(API)
rho_gas_sc = GAS_GR * 1.225
rho_wtr_sc = 1_025.0

denom  = rho_oil_sc + GOR_SC * rho_gas_sc + WOR_SC * rho_wtr_sc
f_oil  = rho_oil_sc          / denom
f_gas  = GOR_SC * rho_gas_sc / denom
f_wtr  = WOR_SC * rho_wtr_sc / denom

q_oil_sc = total_mdot * f_oil / rho_oil_sc * 3600
q_gas_sc = total_mdot * f_gas / rho_gas_sc * 3600
q_wtr_sc = total_mdot * f_wtr / rho_wtr_sc * 3600

print()
print("Surface rates at separator (standard conditions):")
print(f"  Oil:    {q_oil_sc:.1f} m³/h")
print(f"  Gas:    {q_gas_sc:.0f} m³/h")
print(f"  Water:  {q_wtr_sc:.1f} m³/h")
