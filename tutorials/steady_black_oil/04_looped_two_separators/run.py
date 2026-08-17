"""
Tutorial 04 — Black-Oil Looped Network: Two Flow Delivery Points
================================================================
A diamond loop (two parallel paths) with a branch tapped off one of the
junction nodes.  Two delivery points are served simultaneously: the main
separator at the end of the diamond loop, and a satellite separator on
the branch.  The solver distributes total flow between the two parallel
paths of the loop while satisfying both delivery demands.

Geometry
--------
                     PipeA (5 km, D=0.20 m)        PipeC (5 km, D=0.18 m)
      [Node 1] ─────────────────────────> [Node 2] ─────────────────────────> [Node 4]
      P = 8 MPa                                                               ṁ = 80 kg/s
      T = 60 °C                                                               Separator A
          │                                 ↑
          │  PipeB (5 km, D=0.18 m)     PipeD (5 km, D=0.16 m)
          └───────────────────────> [Node 3] ─────────────────────────> Node 4 (loop)
                                       │
                                       │  PipeE (2 km, D=0.14 m)
                                       └──────────────────────────> [Node 5]
                                                                    ṁ = 20 kg/s
                                                                    Separator B

Loop : 1 → 2 → 4 and 1 → 3 → 4  (diamond, two parallel paths)
Branch: Node 3 → Node 5           (flow tap to Separator B)

Total delivery: 80 + 20 = 100 kg/s
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
    # Diamond loop
    Pipe(component_id="pipe_A_1_2", start_node=1, end_node=2,
         diameter_m=0.20, length_m=5_000.0, absolute_roughness_m=46e-6,
         heat_transfer_coefficient_w_per_m2k=5.0, ambient_temperature_c=15.0),
    Pipe(component_id="pipe_B_1_3", start_node=1, end_node=3,
         diameter_m=0.18, length_m=5_000.0, absolute_roughness_m=46e-6,
         heat_transfer_coefficient_w_per_m2k=5.0, ambient_temperature_c=15.0),
    Pipe(component_id="pipe_C_2_4", start_node=2, end_node=4,
         diameter_m=0.18, length_m=5_000.0, absolute_roughness_m=46e-6,
         heat_transfer_coefficient_w_per_m2k=5.0, ambient_temperature_c=15.0),
    Pipe(component_id="pipe_D_3_4", start_node=3, end_node=4,
         diameter_m=0.16, length_m=5_000.0, absolute_roughness_m=46e-6,
         heat_transfer_coefficient_w_per_m2k=5.0, ambient_temperature_c=15.0),
    # Branch to satellite separator
    Pipe(component_id="pipe_E_3_5", start_node=3, end_node=5,
         diameter_m=0.14, length_m=2_000.0, absolute_roughness_m=46e-6,
         heat_transfer_coefficient_w_per_m2k=5.0, ambient_temperature_c=15.0),
]

DELIVERY_4 = 80.0   # kg/s at Node 4 (main separator)
DELIVERY_5 = 20.0   # kg/s at Node 5 (satellite separator)

case = NetworkCase(
    name             = "Black-Oil Looped Network — Two Flow Delivery Points",
    fluid_model      = fluid,
    pressure_inlets  = (PressureBoundary(node_id=1, pressure_pa=8e6),),
    pressure_outlets = (),
    flow_outlets     = (
        FlowBoundary(node_id=4, mass_flow_kg_per_s=DELIVERY_4),
        FlowBoundary(node_id=5, mass_flow_kg_per_s=DELIVERY_5),
    ),
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
for cf in result.component_flows:
    print(f"{cf.label:30s}  {cf.mass_flow_kg_per_s:>10.3f}  {cf.volumetric_flow_m3_per_h:>13.2f}")
print()

# ── Loop split summary ────────────────────────────────────────────────────────
flows = {cf.label.split(":")[-1].strip(): cf.mass_flow_kg_per_s
         for cf in result.component_flows}
upper = flows.get("pipe_A_1_2", 0.0)
lower = flows.get("pipe_B_1_3", 0.0)
branch = flows.get("pipe_E_3_5", 0.0)
total = upper + lower
print(f"Upper path (A→C):  {upper:.3f} kg/s  ({100*upper/total:.1f} % of total)")
print(f"Lower path (B→D/E): {lower:.3f} kg/s  ({100*lower/total:.1f} % of total)")
print(f"  Branch to Sep B (E): {branch:.3f} kg/s")
print(f"Total from inlet:  {total:.3f} kg/s")
print()
print(f"Separator A (Node 4): P = {result.node_pressures_pa[4]/1e6:.3f} MPa")
print(f"Separator B (Node 5): P = {result.node_pressures_pa[5]/1e6:.3f} MPa")
