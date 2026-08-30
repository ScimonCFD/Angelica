"""
Tutorial 05 — Black-Oil Two-Reservoir Blending
===============================================
Two reservoirs with different fluid properties (API gravity, GOR, WOR) feed
a common trunk that delivers to a single separator.  Each inlet carries its
own composition; the solver propagates and mixes them at the junction.

Geometry
--------
  Node 1 ──PipeA── Node 3 ──PipeC── Node 4
  (Res A)  5 km            3 km     (Separator,
  P=9 MPa  D=0.18          D=0.22    P=2 MPa)
  T=70 °C

  Node 2 ──PipeB── Node 3
  (Res B)  5 km
  P=8 MPa  D=0.16
  T=60 °C

Reservoir A — light crude:   32°API, GOR=25 m³/m³, WOR=0.3, gas_gravity=0.65
Reservoir B — heavy crude:   22°API, GOR=10 m³/m³, WOR=1.5, gas_gravity=0.70

The solver mixes both fluids at Node 3 and delivers the blend to Node 4.
The outlet composition is a mass-weighted average of the two inlet streams.
"""
import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[4] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.core.case import InletFluidBC, NetworkCase, PressureBoundary, ThermalBoundary
from angelica.core.components import Pipe
from angelica.properties.black_oil import BlackOilFluid, bubble_point_pa
from angelica.solvers import SteadyBlackOilSolver

# ── Reservoir fluid definitions ───────────────────────────────────────────────
API_A,  GAS_A,  GOR_A,  WOR_A  = 32.0, 0.65, 25.0, 0.3
API_B,  GAS_B,  GOR_B,  WOR_B  = 22.0, 0.70, 10.0, 1.5
T_A, T_B = 70.0, 60.0

Pb_A = bubble_point_pa(GOR_A, GAS_A, API_A, T_A)
Pb_B = bubble_point_pa(GOR_B, GAS_B, API_B, T_B)
print(f"Bubble point  Reservoir A ({API_A}°API, GOR={GOR_A}): {Pb_A/1e6:.2f} MPa")
print(f"Bubble point  Reservoir B ({API_B}°API, GOR={GOR_B}): {Pb_B/1e6:.2f} MPa")
print()

# ── Network ───────────────────────────────────────────────────────────────────
pipes = [
    Pipe(component_id="pipe_A",   start_node=1, end_node=3,
         diameter_m=0.18, length_m=5_000.0, absolute_roughness_m=46e-6,
         heat_transfer_coefficient_w_per_m2k=5.0, ambient_temperature_c=15.0),
    Pipe(component_id="pipe_B",   start_node=2, end_node=3,
         diameter_m=0.16, length_m=5_000.0, absolute_roughness_m=46e-6,
         heat_transfer_coefficient_w_per_m2k=5.0, ambient_temperature_c=15.0),
    Pipe(component_id="pipe_C",   start_node=3, end_node=4,
         diameter_m=0.22, length_m=3_000.0, absolute_roughness_m=46e-6,
         heat_transfer_coefficient_w_per_m2k=5.0, ambient_temperature_c=15.0),
]

# global fluid_model — used as fallback for composition initialisation
default_fluid = BlackOilFluid(
    api_gravity=API_A, gas_gravity=GAS_A,
    gor_sc_m3_per_m3=GOR_A, wor_sc_m3_per_m3=WOR_A,
    reference_pressure_pa=5e6, reference_temperature_c=T_A,
)

case = NetworkCase(
    name             = "Two-Reservoir Blending",
    fluid_model      = default_fluid,
    pressure_inlets  = (
        PressureBoundary(node_id=1, pressure_pa=9e6),
        PressureBoundary(node_id=2, pressure_pa=8e6),
    ),
    pressure_outlets = (PressureBoundary(node_id=4, pressure_pa=2e6),),
    components       = tuple(pipes),
    thermal_inlets   = (
        ThermalBoundary(node_id=1, temperature_c=T_A, bc_type="fixed_temperature"),
        ThermalBoundary(node_id=2, temperature_c=T_B, bc_type="fixed_temperature"),
    ),
    inlet_fluid_bcs  = (
        InletFluidBC(node_id=1, api_gravity=API_A, gas_gravity=GAS_A,
                     gor_sc_m3_per_m3=GOR_A, wor_sc_m3_per_m3=WOR_A),
        InletFluidBC(node_id=2, api_gravity=API_B, gas_gravity=GAS_B,
                     gor_sc_m3_per_m3=GOR_B, wor_sc_m3_per_m3=WOR_B),
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

# ── Node results ──────────────────────────────────────────────────────────────
print(f"{'Node':>6}  {'P (MPa)':>9}  {'T (°C)':>8}")
for node_id in sorted(result.node_pressures_pa):
    P = result.node_pressures_pa[node_id]
    T = result.node_temperatures_c[node_id]
    print(f"{node_id:>6}  {P/1e6:>9.3f}  {T:>8.2f}")
print()

# ── Per-pipe flows ────────────────────────────────────────────────────────────
print(f"{'Pipe':20s}  {'ṁ (kg/s)':>10}  {'Q_mix (m³/h)':>13}")
flows = {}
for cf in result.component_flows:
    name = cf.label.split(":")[-1].strip()
    flows[name] = cf.mass_flow_kg_per_s
    print(f"{name:20s}  {cf.mass_flow_kg_per_s:>10.3f}  {cf.volumetric_flow_m3_per_h:>13.2f}")
print()

# ── Blended composition at outlet ────────────────────────────────────────────
m_A = flows.get("pipe_A", 0.0)
m_B = flows.get("pipe_B", 0.0)
m_tot = m_A + m_B

from angelica.properties.dead_oil import dead_oil_density_kg_per_m3


def surface_rates(m_kg_s, api, gg, gor, wor):
    rho_oil = dead_oil_density_kg_per_m3(api)
    rho_gas = gg * 1.225
    rho_wtr = 1_025.0
    denom   = rho_oil + gor * rho_gas + wor * rho_wtr
    q_oil   = m_kg_s * (rho_oil / denom) / rho_oil * 3600
    q_gas   = m_kg_s * (gor * rho_gas / denom) / rho_gas * 3600
    q_wtr   = m_kg_s * (wor * rho_wtr / denom) / rho_wtr * 3600
    return q_oil, q_gas, q_wtr

q_oil_A, q_gas_A, q_wtr_A = surface_rates(m_A, API_A, GAS_A, GOR_A, WOR_A)
q_oil_B, q_gas_B, q_wtr_B = surface_rates(m_B, API_B, GAS_B, GOR_B, WOR_B)

print("Surface rates at separator (standard conditions):")
print(f"  From Reservoir A ({m_A:.1f} kg/s):")
print(f"    Oil {q_oil_A:.1f} m³/h  |  Gas {q_gas_A:.0f} m³/h  |  Water {q_wtr_A:.1f} m³/h")
print(f"  From Reservoir B ({m_B:.1f} kg/s):")
print(f"    Oil {q_oil_B:.1f} m³/h  |  Gas {q_gas_B:.0f} m³/h  |  Water {q_wtr_B:.1f} m³/h")
print(f"  TOTAL ({m_tot:.1f} kg/s):")
print(f"    Oil {q_oil_A+q_oil_B:.1f} m³/h  |  Gas {q_gas_A+q_gas_B:.0f} m³/h  |  Water {q_wtr_A+q_wtr_B:.1f} m³/h")

# ── Blended API at junction ───────────────────────────────────────────────────
api_blend = (m_A * API_A + m_B * API_B) / m_tot
gor_blend = (m_A * GOR_A + m_B * GOR_B) / m_tot
wor_blend = (m_A * WOR_A + m_B * WOR_B) / m_tot
print()
print("Blended fluid at Node 3 (junction):")
print(f"  API = {api_blend:.1f}°  |  GOR = {gor_blend:.1f} m³/m³  |  WOR = {wor_blend:.2f} m³/m³")
