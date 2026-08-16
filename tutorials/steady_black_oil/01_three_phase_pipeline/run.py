"""
Tutorial: Black-Oil Three-Phase Pipeline
=========================================
A single-pipe crude-oil system with dissolved gas and produced water.

Geometry
--------
  [Node 1] ──── Pipe (10 km, D=0.2 m) ──── [Node 2]
   P = 8 MPa                                 P = 2 MPa
   T = 60 °C

Fluid
-----
  Stock-tank oil : 32 °API, γ_g = 0.65
  GOR            : 25 m³/m³  (at standard conditions)
  WOR            : 0.5 m³/m³ (at standard conditions)

Expected PVT behaviour
----------------------
  Bubble point at 60 °C ≈ 5.5 MPa.
  Inlet (8 MPa > Pb) → undersaturated: all gas dissolved, holdup_gas = 0.
  Outlet (2 MPa < Pb) → two-phase: free gas appears, mixture density drops.

  The outer loop converges by iterating until the mixture density field
  (which changes as gas liberates from solution) stops changing.
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
GOR_SC    = 25.0   # m³/m³ at standard conditions → Pb ≈ 5.5 MPa at 60 °C
WOR_SC    = 0.5    # m³/m³ at standard conditions
T_INLET_C = 60.0

fluid = BlackOilFluid(
    api_gravity             = API,
    gas_gravity             = GAS_GR,
    gor_sc_m3_per_m3        = GOR_SC,
    wor_sc_m3_per_m3        = WOR_SC,
    reference_pressure_pa   = 5e6,
    reference_temperature_c = T_INLET_C,
)

# ── PVT preview ───────────────────────────────────────────────────────────────
Pb = bubble_point_pa(GOR_SC, GAS_GR, API, T_INLET_C)
pvt_inlet  = fluid.pvt(8e6, T_INLET_C)
pvt_outlet = fluid.pvt(2e6, T_INLET_C)

print("=" * 60)
print(f"Bubble point at {T_INLET_C:.0f} °C : {Pb/1e6:.2f} MPa")
print(f"  → Inlet (8 MPa) is {'UNDERSATURATED' if pvt_inlet.undersaturated else 'TWO-PHASE'}")
print(f"  → Outlet (2 MPa) is {'UNDERSATURATED' if pvt_outlet.undersaturated else 'TWO-PHASE'}")
print()
print(f"{'':30s}  {'Inlet (8 MPa)':>14s}  {'Outlet (2 MPa)':>14s}")
print(f"{'Rs (m³/m³)':30s}  {pvt_inlet.rs_m3_per_m3:>14.2f}  {pvt_outlet.rs_m3_per_m3:>14.2f}")
print(f"{'Bo (m³_res/m³_sc)':30s}  {pvt_inlet.bo:>14.4f}  {pvt_outlet.bo:>14.4f}")
print(f"{'z-factor':30s}  {pvt_inlet.z:>14.4f}  {pvt_outlet.z:>14.4f}")
print(f"{'Holdup oil':30s}  {pvt_inlet.holdup_oil:>14.4f}  {pvt_outlet.holdup_oil:>14.4f}")
print(f"{'Holdup gas':30s}  {pvt_inlet.holdup_gas:>14.4f}  {pvt_outlet.holdup_gas:>14.4f}")
print(f"{'Holdup water':30s}  {pvt_inlet.holdup_water:>14.4f}  {pvt_outlet.holdup_water:>14.4f}")
print(f"{'ρ_mixture (kg/m³)':30s}  {pvt_inlet.mixture_density_kg_per_m3:>14.1f}  {pvt_outlet.mixture_density_kg_per_m3:>14.1f}")
print(f"{'μ_mixture (cP)':30s}  {pvt_inlet.mixture_viscosity_pa_s*1e3:>14.3f}  {pvt_outlet.mixture_viscosity_pa_s*1e3:>14.3f}")

# ── Network case ──────────────────────────────────────────────────────────────
pipe = Pipe(
    component_id                        = "pipeline_1_2",
    start_node                          = 1,
    end_node                            = 2,
    diameter_m                          = 0.2,
    length_m                            = 10_000.0,
    absolute_roughness_m                = 46e-6,
    heat_transfer_coefficient_w_per_m2k = 5.0,
    ambient_temperature_c               = 15.0,
)

case = NetworkCase(
    name             = "Three-Phase Black-Oil Pipeline",
    fluid_model      = fluid,
    pressure_inlets  = (PressureBoundary(node_id=1, pressure_pa=8e6),),
    pressure_outlets = (PressureBoundary(node_id=2, pressure_pa=2e6),),
    components       = (pipe,),
    thermal_inlets   = (
        ThermalBoundary(node_id=1, temperature_c=T_INLET_C, bc_type="fixed_temperature"),
    ),
)

# ── Solve ─────────────────────────────────────────────────────────────────────
solver = SteadyBlackOilSolver()
result = solver.solve(case)

print()
print("=" * 60)
print(f"Case: {case.name}")
print(f"Converged:           {result.converged}")
print(f"PVT iterations:      {len(result.density_history)}")

cf = result.component_flows[0]
print()
print(f"Mass flow:   {cf.mass_flow_kg_per_s:.4f} kg/s")
print(f"Vol. flow:   {cf.volumetric_flow_m3_per_h:.3f} m³/h  (mixture at actual conditions)")
print()
print(f"T inlet:     {result.node_temperatures_c[1]:.2f} °C")
print(f"T outlet:    {result.node_temperatures_c[2]:.2f} °C")

# ── Surface rates ─────────────────────────────────────────────────────────────
# The composition at surface is fixed: 1 m³_oil + GOR m³_gas + WOR m³_water.
# Mass fractions are therefore constant and independent of reservoir P, T.
from angelica.properties.dead_oil import dead_oil_density_kg_per_m3

rho_oil_sc = dead_oil_density_kg_per_m3(API)
rho_gas_sc = GAS_GR * 1.225
rho_wtr_sc = 1_025.0

denom  = rho_oil_sc + GOR_SC * rho_gas_sc + WOR_SC * rho_wtr_sc
f_oil  = rho_oil_sc        / denom
f_gas  = GOR_SC * rho_gas_sc / denom
f_wtr  = WOR_SC * rho_wtr_sc / denom

mdot = cf.mass_flow_kg_per_s
q_oil_sc = mdot * f_oil / rho_oil_sc * 3600
q_gas_sc = mdot * f_gas / rho_gas_sc * 3600
q_wtr_sc = mdot * f_wtr / rho_wtr_sc * 3600

print()
print("Surface rates (standard conditions):")
print(f"  Oil:    {q_oil_sc:.2f} m³/h")
print(f"  Gas:    {q_gas_sc:.1f} m³/h  ({q_gas_sc/q_oil_sc:.1f} m³/m³  ← GOR input: {GOR_SC})")
print(f"  Water:  {q_wtr_sc:.2f} m³/h  ({q_wtr_sc/q_oil_sc:.2f} m³/m³  ← WOR input: {WOR_SC})")
