"""
Tutorial 04 — Free Water in a Wet-Gas Pipeline
===============================================
Demonstrates Angelica's immiscible-water model.  When water is included in
the component list the solver detects it automatically, performs the VL flash
on the dry (water-free) normalised composition, and determines how much water
condenses as free liquid using the Wagner / IAPWS-IF97 saturation pressure.

No API change is required: adding "water" to the component list and its mole
fraction to the composition vector is all that is needed.

Geometry
--------
  Node 1 ──── Pipe A (50 km, D = 0.20 m) ──── Node 2 ──── Pipe B (50 km) ──── Node 3
  P = 80 bar, T = 50 °C                                                         P = 20 bar
  (inlet: wet natural gas)

Fluid
-----
  Component        mol%
  --------         ----
  methane          75.0
  ethane            8.0
  propane           4.0
  n-butane          1.0
  water            12.0   ← produced water content typical of a wet-gas well

The pipeline starts warm (50 °C) and cools toward 15 °C ambient.  As the gas
cools and expands, more water condenses from the gas phase.

What this tutorial shows
------------------------
  1. PVT comparison: dry basis vs wet basis flash at the same conditions
  2. Water-phase split at different pressures and temperatures
  3. Full network solve: free_water_fraction per pipe in the result
  4. Mass flow of liquid water estimated from the mole-fraction result
"""
import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[5] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.core.case import InletCompositionBC, NetworkCase, PressureBoundary, ThermalBoundary
from angelica.core.components import Pipe
from angelica.properties.compositional_fluid import CompositionalFluid, _flash_properties
from angelica.properties.free_water import water_psat_pa
from angelica.solvers import SteadyCompositionalSolver

# ── Fluid definition ──────────────────────────────────────────────────────────
COMPONENTS_WET = ["methane", "ethane", "propane", "n-butane", "water"]
ZS_WET         = [0.75, 0.08, 0.04, 0.01, 0.12]   # 12 mol% water

# Dry-basis composition (water removed, HC re-normalised)
ZS_DRY_RAW = ZS_WET[:4]
SUM_DRY    = sum(ZS_DRY_RAW)
COMPONENTS_DRY = COMPONENTS_WET[:4]
ZS_DRY         = [z / SUM_DRY for z in ZS_DRY_RAW]

# ── Section 1: PVT comparison — dry vs wet at inlet conditions ────────────────
P_IN  = 80e5    # Pa  (80 bar)
P_OUT = 20e5    # Pa  (20 bar)
T_IN  = 50.0    # °C  (inlet)
T_AMB = 15.0    # °C  (ambient / outlet approach temperature)

ZS_WET_T = tuple(round(z, 4) for z in ZS_WET)
ZS_DRY_T = tuple(round(z, 4) for z in ZS_DRY)

rho_wet_in, mu_wet_in, _, _, VF_wet_in, fw_in = _flash_properties(
    tuple(COMPONENTS_WET), P_IN, T_IN, ZS_WET_T
)
rho_dry_in, mu_dry_in, _, _, VF_dry_in, _      = _flash_properties(
    tuple(COMPONENTS_DRY), P_IN, T_IN, ZS_DRY_T
)

print("=" * 65)
print("Section 1 — PVT comparison at inlet (80 bar, 50 °C)")
print("-" * 65)
print(f"{'':35s}  {'Wet gas':>12s}  {'Dry gas':>12s}")
print(f"{'ρ  (kg/m³)':35s}  {rho_wet_in:>12.2f}  {rho_dry_in:>12.2f}")
print(f"{'μ  (Pa·s × 10⁻⁵)':35s}  {mu_wet_in*1e5:>12.4f}  {mu_dry_in*1e5:>12.4f}")
print(f"{'VF (vapour fraction)':35s}  {VF_wet_in:>12.4f}  {VF_dry_in:>12.4f}")
print(f"{'free water (mol frac of feed)':35s}  {fw_in:>12.4f}  {'N/A':>12s}")
if fw_in > 0:
    print(f"  → {fw_in*100:.2f}% of the feed condenses as liquid water at inlet")
else:
    print("  → All water remains in vapour phase at inlet conditions")
print()

# ── Section 2: Water-phase split scan across T and P ─────────────────────────
print("=" * 65)
print("Section 2 — Free-water fraction at different T and P")
print("-" * 65)
print(f"{'T (°C)':>8s}  {'P (bar)':>8s}  {'Psat_H2O (bar)':>15s}  {'free water (mol%)':>18s}")
print("-" * 65)

conditions = [
    (50.0,  80.0),
    (30.0,  80.0),
    (15.0,  80.0),
    (50.0,  40.0),
    (30.0,  40.0),
    (15.0,  20.0),
    (15.0,  10.0),
]
for T_c, P_bar in conditions:
    P_pa = P_bar * 1e5
    T_K  = T_c + 273.15
    _, _, _, _, _, fw = _flash_properties(tuple(COMPONENTS_WET), P_pa, T_c, ZS_WET_T)
    Psat = water_psat_pa(T_K)
    print(f"{T_c:>8.1f}  {P_bar:>8.1f}  {Psat/1e5:>15.4f}  {fw*100:>18.3f}")
print()

# ── Section 3: Network solve ──────────────────────────────────────────────────
fluid = CompositionalFluid(components=COMPONENTS_WET, default_zs=ZS_WET)

pipe_A = Pipe(
    component_id         = "pipe_A",
    start_node           = 1,
    end_node             = 2,
    diameter_m           = 0.20,
    length_m             = 50_000.0,
    absolute_roughness_m = 46e-6,
    heat_transfer_coefficient_w_per_m2k = 3.0,
    ambient_temperature_c = T_AMB,
)
pipe_B = Pipe(
    component_id         = "pipe_B",
    start_node           = 2,
    end_node             = 3,
    diameter_m           = 0.20,
    length_m             = 50_000.0,
    absolute_roughness_m = 46e-6,
    heat_transfer_coefficient_w_per_m2k = 3.0,
    ambient_temperature_c = T_AMB,
)

case = NetworkCase(
    name             = "Wet-Gas Pipeline — Free Water",
    fluid_model      = fluid,
    pressure_inlets  = (PressureBoundary(node_id=1, pressure_pa=P_IN),),
    pressure_outlets = (PressureBoundary(node_id=3, pressure_pa=P_OUT),),
    components       = (pipe_A, pipe_B),
    thermal_inlets   = (
        ThermalBoundary(node_id=1, temperature_c=T_IN, bc_type="fixed_temperature"),
    ),
    inlet_composition_bcs = (
        InletCompositionBC(node_id=1, zs=tuple(ZS_WET)),
    ),
)

solver = SteadyCompositionalSolver()
result = solver.solve(case)

print("=" * 65)
print("Section 3 — Network solve results")
print("-" * 65)
print(f"Case:       {case.name}")
print(f"Converged:  {result.converged}")
print()

print(f"{'Node':>5}  {'P (bar)':>9}  {'T (°C)':>8}")
for nid in sorted(result.node_pressures_pa):
    P = result.node_pressures_pa[nid]
    T = result.node_temperatures_c[nid]
    print(f"{nid:>5}  {P/1e5:>9.3f}  {T:>8.2f}")
print()

# ── Section 4: Free-water results per pipe ────────────────────────────────────
print(f"{'Pipe':>6}  {'mdot (kg/s)':>11}  {'VF':>6}  {'free water':>11}  {'water (kg/s est.)':>18}")
print("-" * 65)

# Molecular weight of water (g/mol) and approximate average MW of wet feed
MW_WATER = 18.015
# Rough MW of full wet feed (to convert mol frac → mass frac for estimation)
MW_MIX_APPROX = sum(z * mw for z, mw in zip(
    ZS_WET,
    [16.04, 30.07, 44.10, 58.12, 18.015]  # approximate MWs
))

for cf in result.component_flows:
    mdot = cf.mass_flow_kg_per_s
    vf   = cf.vapor_fraction if cf.vapor_fraction is not None else float("nan")
    fw   = cf.free_water_fraction
    # Estimate liquid water mass flow: n_wl/n_total * (MW_water/MW_mix) * mdot_total
    mf_water = fw * MW_WATER / MW_MIX_APPROX  # approximate mass fraction of liquid water
    mdot_water_est = abs(mdot) * mf_water
    print(
        f"{cf.label:>6}  {mdot:>11.4f}  {vf:>6.3f}  "
        f"{fw*100:>10.3f}%  {mdot_water_est*3600:>14.2f} kg/h"
    )

print()
print("Note: 'free water' is the mole fraction of the feed that is liquid")
print("water at each pipe's average P and T.  The mass flow estimate uses")
print("an approximate mixture molecular weight.")
print()
print(f"Global mass balance error: {result.global_balance.mass_error_pct:.4f} %")
