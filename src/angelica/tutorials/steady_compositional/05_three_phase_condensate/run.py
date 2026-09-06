"""
Tutorial 05 — Three-Phase Rich Condensate Pipeline
===================================================
Demonstrates a three-phase result: Gas + HC liquid + Free water.

At the inlet (80 bar, 110 °C) the rich condensate gas is above its
hydrocarbon dew point, so there are only two phases: gas and free water.
As the fluid cools along the pipeline it crosses the HC dew point
(≈ 105 °C at 80 bar) and a hydrocarbon liquid phase appears.
By the outlet (30 bar, near ambient) the system has three coexisting phases.

Fluid
-----
  Component          mol%
  ─────────          ────
  methane            55.0   lighter carrier gas
  ethane             10.0
  propane            10.0
  n-butane            8.0
  n-pentane           5.0   heavy enough to condense at pipeline T
  n-hexane            5.0
  water               7.0   produced water

Geometry
--------
  Node 1 ── Pipe A (60 km) ── Node 2 ── Pipe B (60 km) ── Node 3
  80 bar, 110 °C                                            30 bar
  (inlet: above HC dew point)     (outlet: three-phase)

What this tutorial shows
------------------------
  1. Phase-state scan: where the HC dew point falls for this mixture
  2. Network solve: V, L, free water, and phase compositions per pipe
  3. Reading the Phase Report from the CSV/Excel export
"""
import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[5] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.core.case import InletCompositionBC, NetworkCase, PressureBoundary, ThermalBoundary
from angelica.core.components import Pipe
from angelica.properties.compositional_fluid import CompositionalFluid, _hc_phase_compositions
from angelica.properties.free_water import find_water_index
from angelica.solvers import SteadyCompositionalSolver

# ── Fluid definition ──────────────────────────────────────────────────────────
COMPONENTS = ["methane", "ethane", "propane", "n-butane", "n-pentane", "n-hexane", "water"]
ZS         = [0.55, 0.10, 0.10, 0.08, 0.05, 0.05, 0.07]

P_IN  = 80e5    # Pa  (80 bar)
P_OUT = 30e5    # Pa  (30 bar)
T_IN  = 110.0   # °C  — above HC dew point at 80 bar
T_AMB = 15.0    # °C  — ambient / soil temperature

# HC-only (water stripped, renormalised)
w_idx     = find_water_index(tuple(COMPONENTS))
HC_NAMES  = tuple(n for i, n in enumerate(COMPONENTS) if i != w_idx)
hc_zs_raw = tuple(z for i, z in enumerate(ZS)         if i != w_idx)
SUM_HC    = sum(hc_zs_raw)
HC_ZS     = tuple(z / SUM_HC for z in hc_zs_raw)

# ── Section 1: HC dew-point scan at 80 bar ───────────────────────────────────
print("=" * 68)
print("Section 1 — HC phase state scan at 80 bar (where does liquid appear?)")
print("-" * 68)
print(f"{'T (°C)':>7}  {'VF_hc':>8}  {'V  gas':>8}  {'L  HCliq':>10}  {'Phase (HC)':}")
print("-" * 68)

for T_c in [120.0, 115.0, 110.0, 105.0, 100.0, 90.0, 70.0, 50.0, 30.0]:
    P_pa = P_IN
    _, _, VF_hc = _hc_phase_compositions(HC_NAMES, P_pa, T_c, HC_ZS)
    V    = VF_hc * SUM_HC
    L_hc = (1.0 - VF_hc) * SUM_HC
    phase = "Single-phase gas" if L_hc < 1e-4 else f"Gas + HC liquid  (L = {L_hc:.4f})"
    print(f"{T_c:>7.1f}  {VF_hc:>8.4f}  {V:>8.4f}  {L_hc:>10.4f}  {phase}")

print()

# ── Section 2: Network solve ──────────────────────────────────────────────────
fluid = CompositionalFluid(components=COMPONENTS, default_zs=ZS)

pipe_A = Pipe(
    component_id                        = "pipe_A",
    start_node                          = 1,
    end_node                            = 2,
    diameter_m                          = 0.25,
    length_m                            = 60_000.0,
    absolute_roughness_m                = 46e-6,
    heat_transfer_coefficient_w_per_m2k = 3.0,
    ambient_temperature_c               = T_AMB,
)
pipe_B = Pipe(
    component_id                        = "pipe_B",
    start_node                          = 2,
    end_node                            = 3,
    diameter_m                          = 0.25,
    length_m                            = 60_000.0,
    absolute_roughness_m                = 46e-6,
    heat_transfer_coefficient_w_per_m2k = 3.0,
    ambient_temperature_c               = T_AMB,
)

case = NetworkCase(
    name              = "Three-Phase Rich Condensate Pipeline",
    fluid_model       = fluid,
    pressure_inlets   = (PressureBoundary(node_id=1, pressure_pa=P_IN),),
    pressure_outlets  = (PressureBoundary(node_id=3, pressure_pa=P_OUT),),
    components        = (pipe_A, pipe_B),
    thermal_inlets    = (ThermalBoundary(node_id=1, temperature_c=T_IN, bc_type="fixed_temperature"),),
    inlet_composition_bcs = (InletCompositionBC(node_id=1, zs=tuple(ZS)),),
)

solver = SteadyCompositionalSolver()
result = solver.solve(case)

print("=" * 68)
print("Section 2 — Network solve results")
print("-" * 68)
print(f"Case:       {case.name}")
print(f"Converged:  {result.converged}")
print()

print(f"{'Node':>5}  {'P (bar)':>9}  {'T (°C)':>8}")
for nid in sorted(result.node_pressures_pa):
    P = result.node_pressures_pa[nid]
    T = result.node_temperatures_c.get(nid, float("nan"))
    print(f"{nid:>5}  {P/1e5:>9.3f}  {T:>8.2f}")
print()

# Phase report per pipe
print(f"{'Pipe':12}  {'T_in':>7}  {'T_out':>7}  {'V gas':>8}  {'L HCliq':>9}  {'FW':>8}  {'Phase state'}")
print("-" * 75)
for cf in result.component_flows:
    vf = cf.vapor_fraction if cf.vapor_fraction is not None else float("nan")
    lf = cf.liquid_fraction if cf.liquid_fraction is not None else 0.0
    fw = cf.free_water_fraction
    phases = []
    if vf  > 1e-3: phases.append("Gas")
    if lf  > 1e-3: phases.append("HC liq")
    if fw  > 1e-3: phases.append("Free water")
    label = " + ".join(phases) if phases else "—"
    t_in  = f"{cf.temperature_in_c:>7.2f}"  if cf.temperature_in_c  is not None else "    —  "
    t_out = f"{cf.temperature_out_c:>7.2f}" if cf.temperature_out_c is not None else "    —  "
    print(f"{cf.label:12}  {t_in}  {t_out}  {vf:>8.4f}  {lf:>9.4f}  {fw:>8.4f}  {label}")

print()

# Phase compositions
print("Gas phase composition  y_i  (within gas phase):")
hc_names_display = [n for n in COMPONENTS if n != "water"]
print(f"  {'Pipe':12}  " + "  ".join(f"{n:>10}" for n in hc_names_display))
for cf in result.component_flows:
    y = cf.gas_phase_zs
    if y:
        vals = "  ".join(f"{v:>10.4f}" for v in y)
        print(f"  {cf.label:12}  {vals}")

print()
has_liq = any((cf.liquid_fraction or 0.0) > 1e-3 for cf in result.component_flows)
if has_liq:
    print("HC liquid phase composition  x_i  (within HC liquid phase):")
    print(f"  {'Pipe':12}  " + "  ".join(f"{n:>10}" for n in hc_names_display))
    for cf in result.component_flows:
        x = cf.liquid_phase_zs
        if x and (cf.liquid_fraction or 0.0) > 1e-3:
            vals = "  ".join(f"{v:>10.4f}" for v in x)
            print(f"  {cf.label:12}  {vals}")
else:
    print("No HC liquid phase at these conditions.")

print()
print(f"Global mass balance error: {result.global_balance.mass_error_pct:.4f} %")
