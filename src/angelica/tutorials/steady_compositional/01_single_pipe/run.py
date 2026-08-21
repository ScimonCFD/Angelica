"""
Tutorial 01 — Compositional Single Pipe
========================================
A single gas pipe transporting a methane/ethane binary mixture.  The EOS
flash is performed at each iteration using the ``thermo`` library, showing
how density and viscosity change along the pressure gradient.

Geometry
--------
  Node 1 ──── Pipe (10 km, D = 0.15 m) ──── Node 2
  P = 100 bar                                P = 20 bar
  T = 60 °C  (fixed inlet temperature)

Fluid
-----
  Components :  methane (CH₄) 80 mol%,  ethane (C₂H₆) 20 mol%
  Conditions :  well above the dew point at both ends → single-phase gas

What the outer loop does
------------------------
  Each outer iteration:
    1. Flash at (P, T, z) for every pipe → density, viscosity, Cp, k
    2. SIMPLE hydraulic solve (inner loop)
    3. Energy equation → update temperatures
    4. Converge when max|Δρ/ρ| and max|ΔT| fall below tolerance
"""
import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[4] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.core.case import NetworkCase, PressureBoundary, ThermalBoundary, InletCompositionBC
from angelica.core.components import Pipe
from angelica.properties.compositional_fluid import CompositionalFluid, _flash_properties
from angelica.solvers import SteadyCompositionalSolver

# ── Fluid definition ──────────────────────────────────────────────────────────
COMPONENTS = ["methane", "ethane"]
ZS         = [0.80, 0.20]   # overall mole fractions
T_INLET_C  = 60.0

fluid = CompositionalFluid(components=COMPONENTS, default_zs=ZS)

# ── PVT preview: flash at inlet and outlet conditions ─────────────────────────
P_IN  = 100e5   # Pa
P_OUT =  20e5   # Pa
ZS_T  = tuple(round(z, 4) for z in ZS)

rho_in,  mu_in,  Cp_in,  k_in  = _flash_properties(tuple(COMPONENTS), P_IN,  T_INLET_C, ZS_T)
rho_out, mu_out, Cp_out, k_out = _flash_properties(tuple(COMPONENTS), P_OUT, T_INLET_C, ZS_T)

print("=" * 60)
print("EOS flash preview (thermo library)")
print("-" * 60)
print(f"{'':30s}  {'Inlet (100 bar)':>14s}  {'Outlet (20 bar)':>14s}")
print(f"{'ρ  (kg/m³)':30s}  {rho_in:>14.2f}  {rho_out:>14.2f}")
print(f"{'μ  (Pa·s × 10⁻⁵)':30s}  {mu_in*1e5:>14.4f}  {mu_out*1e5:>14.4f}")
print(f"{'Cp (J/(kg·K))':30s}  {Cp_in:>14.1f}  {Cp_out:>14.1f}")
print(f"{'k  (W/(m·K))':30s}  {k_in*1e3:>14.4f}  {k_out*1e3:>14.4f}")
print()

# ── Network case ──────────────────────────────────────────────────────────────
pipe = Pipe(
    component_id         = "gas_pipe",
    start_node           = 1,
    end_node             = 2,
    diameter_m           = 0.15,
    length_m             = 10_000.0,
    absolute_roughness_m = 46e-6,           # commercial steel
    heat_transfer_coefficient_w_per_m2k = 5.0,
    ambient_temperature_c = 15.0,
)

case = NetworkCase(
    name             = "Compositional Single Pipe",
    fluid_model      = fluid,
    pressure_inlets  = (PressureBoundary(node_id=1, pressure_pa=P_IN),),
    pressure_outlets = (PressureBoundary(node_id=2, pressure_pa=P_OUT),),
    components       = (pipe,),
    thermal_inlets   = (
        ThermalBoundary(node_id=1, temperature_c=T_INLET_C, bc_type="fixed_temperature"),
    ),
    inlet_composition_bcs = (
        InletCompositionBC(node_id=1, zs=tuple(ZS)),
    ),
)

# ── Solve ─────────────────────────────────────────────────────────────────────
solver = SteadyCompositionalSolver()
result = solver.solve(case)

cf = result.component_flows[0]

print("=" * 60)
print(f"Case:              {case.name}")
print(f"Converged:         {result.converged}")
print(f"Outer iterations:  {len(result.density_history)}")
print()

print(f"{'Node':>6}  {'P (bar)':>9}  {'T (°C)':>8}")
for nid in sorted(result.node_pressures_pa):
    P = result.node_pressures_pa[nid]
    T = result.node_temperatures_c[nid]
    print(f"{nid:>6}  {P/1e5:>9.3f}  {T:>8.2f}")
print()

print(f"Mass flow:       {cf.mass_flow_kg_per_s:.4f} kg/s")
print(f"Volumetric flow: {cf.volumetric_flow_m3_per_h:.2f} m³/h  (at pipe conditions)")
print()

# Approximate velocity at inlet using inlet density
v_inlet = cf.mass_flow_kg_per_s / (rho_in * pipe.area_m2)
v_outlet = cf.mass_flow_kg_per_s / (rho_out * pipe.area_m2)
print(f"Velocity @ inlet:  {v_inlet:.2f} m/s")
print(f"Velocity @ outlet: {v_outlet:.2f} m/s  (gas expands as P drops)")
print()

print(f"Global mass balance error: {result.global_balance.mass_error_pct:.4f} %")
