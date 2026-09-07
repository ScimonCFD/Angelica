"""
Tutorial 06 — Three-Phase Rich Condensate Loop
===============================================
Same rich condensate fluid as tutorial 05 (Gas + HC liquid + Free water),
now in a looped (ring) network.  Flow leaves the source, reaches a split
junction, travels through two independent paths of different resistance,
recombines at a merge junction, and exits at the sink.

The loop forces the solver to find a flow distribution that satisfies both
mass balance at every node AND the same pressure at the merge junction
regardless of which path is taken.

Geometry
--------

  Source (Node 1, 80 bar, 110 °C)
         |
     Pipe A  (20 km, D = 0.25 m)
         |
     Split junction (Node 2)
       /                     \\
  Pipe B (upper)          Pipe C (lower-left)
  60 km, D = 0.20 m       30 km, D = 0.25 m
  (long / narrow)              |
       \\               Lower junction (Node 3)
        \\                     |
         \\              Pipe D (lower-right)
          \\             30 km, D = 0.25 m
           \\                  |
     Merge junction (Node 4)
         |
     Pipe E  (20 km, D = 0.25 m)
         |
     Sink (Node 5, 30 bar)

Upper path: Pipe B alone  (60 km, D = 0.20)  — longer and narrower → more resistance
Lower path: Pipe C + Pipe D  (60 km total, D = 0.25)  — larger diameter → less resistance
→ more mass flow takes the lower path

Fluid
-----
  Component          mol%
  ─────────          ────
  methane            55.0
  ethane             10.0
  propane            10.0
  n-butane            8.0
  n-pentane           5.0
  n-hexane            5.0
  water               7.0

What this tutorial shows
------------------------
  1. Hydraulic loop: the SIMPLE algorithm finds how flow distributes
     between the upper path (Pipe B: long, narrow → more resistance) and
     the lower path (Pipes C+D: same total length but larger diameter →
     less resistance).  More flow takes the lower path.
  2. Different phase fractions on each path: Pipe B runs cooler (longer,
     more heat loss) so it accumulates more HC liquid condensate.
  3. Gas/liquid compositions (y_i, x_i) in each pipe, demonstrating
     K-value separation between light (methane) and heavy (hexane)
     components.
"""
import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[5] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.core.case import InletCompositionBC, NetworkCase, PressureBoundary, ThermalBoundary
from angelica.core.components import Pipe
from angelica.properties.compositional_fluid import CompositionalFluid
from angelica.solvers import SteadyCompositionalSolver

# ── Fluid ─────────────────────────────────────────────────────────────────────
COMPONENTS = ["methane", "ethane", "propane", "n-butane", "n-pentane", "n-hexane", "water"]
ZS         = [0.55, 0.10, 0.10, 0.08, 0.05, 0.05, 0.07]

P_IN  = 80e5    # Pa  (80 bar)
P_OUT = 30e5    # Pa  (30 bar)
T_IN  = 110.0   # °C
T_AMB = 15.0    # °C

fluid = CompositionalFluid(components=COMPONENTS, default_zs=ZS)

# ── Network ───────────────────────────────────────────────────────────────────
#   1 → [A] → 2 → [B] ──────────→ 4 → [E] → 5
#              ↘ [C] → 3 → [D] ↗
pipe_A = Pipe(
    component_id="pipe_A", start_node=1, end_node=2,
    diameter_m=0.25, length_m=20_000.0, absolute_roughness_m=46e-6,
    heat_transfer_coefficient_w_per_m2k=3.0, ambient_temperature_c=T_AMB,
)
pipe_B = Pipe(
    component_id="pipe_B", start_node=2, end_node=4,
    diameter_m=0.20, length_m=60_000.0, absolute_roughness_m=46e-6,
    heat_transfer_coefficient_w_per_m2k=3.0, ambient_temperature_c=T_AMB,
)
pipe_C = Pipe(
    component_id="pipe_C", start_node=2, end_node=3,
    diameter_m=0.25, length_m=30_000.0, absolute_roughness_m=46e-6,
    heat_transfer_coefficient_w_per_m2k=3.0, ambient_temperature_c=T_AMB,
)
pipe_D = Pipe(
    component_id="pipe_D", start_node=3, end_node=4,
    diameter_m=0.25, length_m=30_000.0, absolute_roughness_m=46e-6,
    heat_transfer_coefficient_w_per_m2k=3.0, ambient_temperature_c=T_AMB,
)
pipe_E = Pipe(
    component_id="pipe_E", start_node=4, end_node=5,
    diameter_m=0.25, length_m=20_000.0, absolute_roughness_m=46e-6,
    heat_transfer_coefficient_w_per_m2k=3.0, ambient_temperature_c=T_AMB,
)

case = NetworkCase(
    name              = "Three-Phase Rich Condensate Loop",
    fluid_model       = fluid,
    pressure_inlets   = (PressureBoundary(node_id=1, pressure_pa=P_IN),),
    pressure_outlets  = (PressureBoundary(node_id=5, pressure_pa=P_OUT),),
    components        = (pipe_A, pipe_B, pipe_C, pipe_D, pipe_E),
    thermal_inlets    = (ThermalBoundary(node_id=1, temperature_c=T_IN, bc_type="fixed_temperature"),),
    inlet_composition_bcs = (InletCompositionBC(node_id=1, zs=tuple(ZS)),),
)

# ── Solve ─────────────────────────────────────────────────────────────────────
solver = SteadyCompositionalSolver()
result = solver.solve(case)

print("=" * 70)
print(f"Case:       {case.name}")
print(f"Converged:  {result.converged}")
print()

# Node pressures and temperatures
labels = {1: "Source", 2: "Split junction", 3: "Lower junction", 4: "Merge junction", 5: "Sink"}
print(f"{'Node':>5}  {'Label':20}  {'P (bar)':>9}  {'T (°C)':>8}")
print("-" * 50)
for nid in sorted(result.node_pressures_pa):
    P = result.node_pressures_pa[nid]
    T = result.node_temperatures_c.get(nid, float("nan"))
    print(f"{nid:>5}  {labels.get(nid,''):20}  {P/1e5:>9.3f}  {T:>8.2f}")
print()

# Flow and phase summary per pipe
print(f"{'Pipe':8}  {'ṁ (kg/s)':>9}  {'T_in':>7}  {'T_out':>7}  "
      f"{'V gas':>7}  {'L HCliq':>9}  {'Free water':>11}  Phase state")
print("-" * 88)
for cf in result.component_flows:
    vf = cf.vapor_fraction  if cf.vapor_fraction  is not None else float("nan")
    lf = cf.liquid_fraction if cf.liquid_fraction is not None else 0.0
    fw = cf.free_water_fraction
    phases = []
    if vf  > 1e-3: phases.append("Gas")
    if lf  > 1e-3: phases.append("HC liq")
    if fw  > 1e-3: phases.append("Free water")
    state = " + ".join(phases) if phases else "—"
    tin  = f"{cf.temperature_in_c:>7.2f}"  if cf.temperature_in_c  is not None else "    —  "
    tout = f"{cf.temperature_out_c:>7.2f}" if cf.temperature_out_c is not None else "    —  "
    name = cf.label.split("_")[-1] if "_" in cf.label else cf.label
    print(f"{name:8}  {cf.mass_flow_kg_per_s:>9.4f}  {tin}  {tout}  "
          f"{vf:>7.4f}  {lf:>9.4f}  {fw:>11.4f}  {state}")

print()
print("Note: Upper path (Pipe B, D=0.20 m) carries less flow than lower path (Pipes C+D, D=0.25 m).")
print("      Pipe B runs cooler (longer) → more HC liquid condensate (L) than C or D.")
print()

# Gas and liquid compositions
hc_names = [c for c in COMPONENTS if c != "water"]
print("Gas phase composition  y_i  (mole fractions within gas phase):")
print(f"  {'Pipe':8}  " + "  ".join(f"{n:>10}" for n in hc_names))
for cf in result.component_flows:
    y = cf.gas_phase_zs
    if y:
        name = cf.label.split("_")[-1] if "_" in cf.label else cf.label
        print(f"  {name:8}  " + "  ".join(f"{v:>10.4f}" for v in y))

print()
has_liq = any((cf.liquid_fraction or 0.0) > 1e-3 for cf in result.component_flows)
if has_liq:
    print("HC liquid phase composition  x_i  (mole fractions within HC liquid phase):")
    print(f"  {'Pipe':8}  " + "  ".join(f"{n:>10}" for n in hc_names))
    for cf in result.component_flows:
        if (cf.liquid_fraction or 0.0) > 1e-3:
            x = cf.liquid_phase_zs
            name = cf.label.split("_")[-1] if "_" in cf.label else cf.label
            if x:
                print(f"  {name:8}  " + "  ".join(f"{v:>10.4f}" for v in x))

print()
print(f"Global mass balance error: {result.global_balance.mass_error_pct:.4f} %")
