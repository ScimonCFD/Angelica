"""
Tutorial 03 — Black-Oil Looped Network with Flow Outlet
========================================================
Trunk-fed diamond loop where the outlet specifies a **mass-flow demand**
instead of a pressure.  The solver finds the outlet pressure that satisfies
the 100 kg/s demand.

Geometry (identical to tutorial 02 — only the outlet BC changes)
--------
                        PipeB  5 km, D=0.20 m             PipeC  5 km, D=0.18 m
                  ┌───────────────────────── Node 3 ─────────────────────────┐
                  │                                                           │
  Node 1 ──PipeA─ Node 2                                                 Node 5 ──PipeF── Node 6
  P=8 MPa  2 km   T-split                                                T-merge   2 km   ṁ=100 kg/s
  T=60 °C  0.22m  │                                                           │    0.22m  (flow BC)
                  └───────────────────────── Node 4 ─────────────────────────┘
                        PipeD  8 km, D=0.15 m             PipeE  3 km, D=0.20 m

  Nodes 2, 5 : T-junctions (bifurcation / convergence)
  Nodes 3, 4 : upper and lower loop junctions
  Node 1 : P = 8 MPa, T = 60 °C (pressure inlet)
  Node 6 : ṁ = 100 kg/s          (flow outlet — solver finds outlet pressure)
"""
import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[4] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.core.case import FlowBoundary, NetworkCase, PressureBoundary, ThermalBoundary
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
    Pipe(component_id="pipe_A_trunk_in",   start_node=1, end_node=2,
         diameter_m=0.22, length_m=2_000.0, absolute_roughness_m=46e-6,
         heat_transfer_coefficient_w_per_m2k=5.0, ambient_temperature_c=15.0),
    Pipe(component_id="pipe_B_upper_L",    start_node=2, end_node=3,
         diameter_m=0.20, length_m=5_000.0, absolute_roughness_m=46e-6,
         heat_transfer_coefficient_w_per_m2k=5.0, ambient_temperature_c=15.0),
    Pipe(component_id="pipe_C_upper_R",    start_node=3, end_node=5,
         diameter_m=0.18, length_m=5_000.0, absolute_roughness_m=46e-6,
         heat_transfer_coefficient_w_per_m2k=5.0, ambient_temperature_c=15.0),
    Pipe(component_id="pipe_D_lower_L",    start_node=2, end_node=4,
         diameter_m=0.15, length_m=8_000.0, absolute_roughness_m=46e-6,
         heat_transfer_coefficient_w_per_m2k=5.0, ambient_temperature_c=15.0),
    Pipe(component_id="pipe_E_lower_R",    start_node=4, end_node=5,
         diameter_m=0.20, length_m=3_000.0, absolute_roughness_m=46e-6,
         heat_transfer_coefficient_w_per_m2k=5.0, ambient_temperature_c=15.0),
    Pipe(component_id="pipe_F_trunk_out",  start_node=5, end_node=6,
         diameter_m=0.22, length_m=2_000.0, absolute_roughness_m=46e-6,
         heat_transfer_coefficient_w_per_m2k=5.0, ambient_temperature_c=15.0),
]

case = NetworkCase(
    name             = "Black-Oil Loop — Flow Outlet",
    fluid_model      = fluid,
    pressure_inlets  = (PressureBoundary(node_id=1, pressure_pa=8e6),),
    pressure_outlets = (),
    flow_outlets     = (FlowBoundary(node_id=6, mass_flow_kg_per_s=100.0),),
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
print(f"{'Pipe':35s}  {'ṁ (kg/s)':>10}  {'Q_mix (m³/h)':>13}")
for cf in result.component_flows:
    print(f"{cf.label:35s}  {cf.mass_flow_kg_per_s:>10.3f}  {cf.volumetric_flow_m3_per_h:>13.2f}")
print()

# ── Loop split ────────────────────────────────────────────────────────────────
flows = {cf.label.split(":")[-1].strip(): cf.mass_flow_kg_per_s
         for cf in result.component_flows}
upper = flows.get("pipe_B_upper_L", 0.0)
lower = flows.get("pipe_D_lower_L", 0.0)
total = upper + lower
P_out = result.node_pressures_pa[6]
print(f"Outlet pressure found: {P_out/1e6:.3f} MPa")
print(f"Upper loop (B+C): {upper:.3f} kg/s  ({100*upper/total:.1f} %)")
print(f"Lower loop (D+E): {lower:.3f} kg/s  ({100*lower/total:.1f} %)")
print(f"Total throughput: {total:.3f} kg/s")
