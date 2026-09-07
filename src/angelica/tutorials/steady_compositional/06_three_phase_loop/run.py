"""
Tutorial 06 — Two-Source Compositional Loop
============================================
Two gas streams with different compositions enter a shared mixing junction.
The blended fluid then travels through a hydraulic loop before reaching the
outlet.  The tutorial shows how the solver handles composition tracking
through mixing, and how the phase behaviour of the combined stream differs
from either source stream on its own.

Source A  (rich condensate, 80 bar, 110 °C)
  55 % CH₄, 10 % C₂H₆, 10 % C₃H₈, 8 % nC₄, 5 % nC₅, 5 % nC₆, 7 % H₂O
  → three-phase in Pipe A  (gas + HC liquid + free water)

Source B  (lean dry gas, 75 bar, 50 °C)
  90 % CH₄, 6 % C₂H₆, 3 % C₃H₈, 1 % nC₄  (no water, no pentane/hexane)
  → single-phase gas in Pipe F

After mixing, the blended composition is intermediate between the two sources:
more methane than the rich condensate alone, less than the lean gas alone.
The water content is diluted but still sufficient for free-water condensation,
so the combined stream remains three-phase through the entire loop.

Geometry
--------

  Source A (Node 1, 80 bar, 110 °C)  ──[Pipe A, 20 km, D=0.25]──┐
                                                                   ├── Mixer (Node 2)
  Source B (Node 7, 75 bar,  50 °C)  ──[Pipe F, 15 km, D=0.20]──┘
                                                                        │
                                                                   [Pipe B, 20 km, D=0.25]
                                                                        │
                                                               Split junction (Node 3)
                                                              /                         \\
                                                [Pipe C, 60 km, D=0.20]     [Pipe D, 30 km, D=0.25]
                                                (upper: long/narrow)          (lower-left)
                                                              \\                         /
                                                               \\              Lower (Node 4)
                                                                \\                  /
                                                                 \\    [Pipe E, 30 km, D=0.25]
                                                                  \\              /
                                                               Merge junction (Node 5)
                                                                        │
                                                                   [Pipe G, 20 km, D=0.25]
                                                                        │
                                                                   Sink (Node 6, 30 bar)

What this tutorial shows
------------------------
  1. Composition mixing: two streams of different composition enter a junction.
     The solver computes the molar-flow-weighted blend automatically.
  2. Phase behaviour contrast: Source A is 3-phase; Source B is single-phase
     dry gas.  After mixing, all downstream pipes are 3-phase because the
     water from Source A persists in the blend.
  3. Hydraulic loop: the SIMPLE algorithm distributes the blended flow between
     the upper path (Pipe C: long, narrow → more resistance) and the lower
     path (Pipes D+E: same total length, larger diameter → less resistance).
  4. K-value separation: y_i (gas) is methane-rich; x_i (HC liquid) is
     enriched in heavier components.
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

ZS_RICH = [0.55, 0.10, 0.10, 0.08, 0.05, 0.05, 0.07]   # Source A — rich condensate
ZS_LEAN = [0.90, 0.06, 0.03, 0.01, 0.00, 0.00, 0.00]   # Source B — lean dry gas
ZS_DEF  = [0.70, 0.08, 0.07, 0.05, 0.03, 0.03, 0.04]   # initial guess for pipes (rough blend)

fluid = CompositionalFluid(components=COMPONENTS, default_zs=ZS_DEF)

T_AMB = 15.0   # °C ambient temperature

# ── Network ───────────────────────────────────────────────────────────────────
def _pipe(cid, s, e, d, L):
    return Pipe(
        component_id=cid, start_node=s, end_node=e,
        diameter_m=d, length_m=L, absolute_roughness_m=46e-6,
        heat_transfer_coefficient_w_per_m2k=3.0, ambient_temperature_c=T_AMB,
    )

pipe_A = _pipe("pipe_A", 1, 2, 0.25, 20_000)   # Source A → Mixer
pipe_F = _pipe("pipe_F", 7, 2, 0.20, 15_000)   # Source B → Mixer
pipe_B = _pipe("pipe_B", 2, 3, 0.25, 20_000)   # Mixer    → Split
pipe_C = _pipe("pipe_C", 3, 5, 0.20, 60_000)   # upper loop path (long, narrow)
pipe_D = _pipe("pipe_D", 3, 4, 0.25, 30_000)   # lower-left
pipe_E = _pipe("pipe_E", 4, 5, 0.25, 30_000)   # lower-right
pipe_G = _pipe("pipe_G", 5, 6, 0.25, 20_000)   # Merge    → Sink

case = NetworkCase(
    name              = "Two-Source Compositional Loop",
    fluid_model       = fluid,
    pressure_inlets   = (
        PressureBoundary(node_id=1, pressure_pa=80e5),   # Source A
        PressureBoundary(node_id=7, pressure_pa=75e5),   # Source B
    ),
    pressure_outlets  = (PressureBoundary(node_id=6, pressure_pa=30e5),),
    components        = (pipe_A, pipe_F, pipe_B, pipe_C, pipe_D, pipe_E, pipe_G),
    thermal_inlets    = (
        ThermalBoundary(node_id=1, temperature_c=110.0, bc_type="fixed_temperature"),
        ThermalBoundary(node_id=7, temperature_c= 50.0, bc_type="fixed_temperature"),
    ),
    inlet_composition_bcs = (
        InletCompositionBC(node_id=1, zs=tuple(ZS_RICH)),
        InletCompositionBC(node_id=7, zs=tuple(ZS_LEAN)),
    ),
)

# ── Solve ─────────────────────────────────────────────────────────────────────
solver = SteadyCompositionalSolver()
result = solver.solve(case)

print("=" * 70)
print(f"Case:       {case.name}")
print(f"Converged:  {result.converged}")
print()

# Node pressures and temperatures
NODE_LABELS = {
    1: "Source A (rich)",
    7: "Source B (lean)",
    2: "Mixer junction",
    3: "Split junction",
    4: "Lower junction",
    5: "Merge junction",
    6: "Sink",
}
print(f"{'Node':>5}  {'Label':22}  {'P (bar)':>9}  {'T (°C)':>8}")
print("-" * 54)
for nid in sorted(result.node_pressures_pa):
    P = result.node_pressures_pa[nid]
    T = result.node_temperatures_c.get(nid, float("nan"))
    print(f"{nid:>5}  {NODE_LABELS.get(nid,''):22}  {P/1e5:>9.3f}  {T:>8.2f}")
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
print("  Source A (rich condensate): 3-phase gas entering Pipe A")
print("  Source B (lean dry gas):    single-phase gas entering Pipe F")
print("  After mixing (Pipe B onward): the blended stream is 3-phase")
print("  (water from Source A survives dilution and still condenses)")
print()

# Gas and liquid compositions
HC_NAMES = [c for c in COMPONENTS if c != "water"]
print("Gas phase composition  y_i  (mole fractions within gas phase):")
print(f"  {'Pipe':8}  " + "  ".join(f"{n:>10}" for n in HC_NAMES))
for cf in result.component_flows:
    y = cf.gas_phase_zs
    if y:
        name = cf.label.split("_")[-1] if "_" in cf.label else cf.label
        print(f"  {name:8}  " + "  ".join(f"{v:>10.4f}" for v in y))

print()
has_liq = any((cf.liquid_fraction or 0.0) > 1e-3 for cf in result.component_flows)
if has_liq:
    print("HC liquid phase composition  x_i  (mole fractions within HC liquid phase):")
    print(f"  {'Pipe':8}  " + "  ".join(f"{n:>10}" for n in HC_NAMES))
    for cf in result.component_flows:
        if (cf.liquid_fraction or 0.0) > 1e-3:
            x = cf.liquid_phase_zs
            name = cf.label.split("_")[-1] if "_" in cf.label else cf.label
            if x:
                print(f"  {name:8}  " + "  ".join(f"{v:>10.4f}" for v in x))

print()
print(f"Global mass balance error: {result.global_balance.mass_error_pct:.4f} %")
