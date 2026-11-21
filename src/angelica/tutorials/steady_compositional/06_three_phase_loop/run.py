"""
Tutorial 06 — Compositional Ring Main: Who Gets the Condensate?
================================================================
Two gas fields feed a ring main that supplies two customers.  Each source
enters the ring at one corner and each customer draws from another corner,
so every customer is fed by both sources along two different paths.  The
hydraulics alone decide how much of each source reaches each customer — and
therefore what composition, how much hydrocarbon liquid and how much free
water each one receives.

Source A  (rich condensate, 80 bar, 110 °C)
  55 % CH₄, 10 % C₂H₆, 10 % C₃H₈, 8 % nC₄, 5 % nC₅, 5 % nC₆, 7 % H₂O
  → gas at the wellhead (above the cricondentherm); it cools in the
    feed line and becomes three-phase (gas + HC liquid + free water)

Source B  (lean dry gas, 75 bar, 50 °C)
  90 % CH₄, 6 % C₂H₆, 3 % C₃H₈, 1 % nC₄  (no water, no C₅₊)
  → single-phase gas everywhere

Geometry
--------

  Source A (1) ──Feed A, 20 km──► N2 ════ Ring north, 25 km, D=0.25 ════► N3 ──5 km──► Plant (6)
  80 bar, 110 °C                   ║                                     ▲            30 bar
                          Ring west, 25 km                      Ring east, 25 km
                              D=0.20                                 D=0.20
                                   ▼                                     ║
  Power station (8) ◄──5 km── N5 ◄════ Ring south, 25 km, D=0.25 ════ N4 ◄──Feed B, 15 km── Source B (7)
  40 bar                                                                                    75 bar, 50 °C

  Source A enters at N2, Source B at N4 (opposite corners).
  The gas plant (30 bar) draws at N3, the power station (40 bar) at N5.

What this tutorial shows
------------------------
  1. A meshed network: each customer is reached from both sources, and each
     source supplies both customers.  The split is not prescribed — it comes
     out of the SIMPLE pressure–flow solution.
  2. Composition mixing at two junctions (N3 and N5) with different blend
     ratios, so the two customers receive different gases.
  3. Phase behaviour: Source A is single-phase at the wellhead but drops into
     the two-phase envelope as it cools; the lean gas never does.  The blend
     each customer receives carries a different share of HC liquid and water.
"""
import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[4]
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

pipes = (
    _pipe("feed_A",     1, 2, 0.25, 20_000),   # Source A → ring (N2)
    _pipe("feed_B",     7, 4, 0.25, 15_000),   # Source B → ring (N4)
    _pipe("ring_north", 2, 3, 0.25, 25_000),
    _pipe("ring_east",  4, 3, 0.20, 25_000),
    _pipe("ring_west",  2, 5, 0.20, 25_000),
    _pipe("ring_south", 4, 5, 0.25, 25_000),
    _pipe("to_plant",   3, 6, 0.25,  5_000),   # N3 → gas plant
    _pipe("to_power",   5, 8, 0.20,  5_000),   # N5 → power station
)

case = NetworkCase(
    name              = "Compositional Ring Main",
    fluid_model       = fluid,
    pressure_inlets   = (
        PressureBoundary(node_id=1, pressure_pa=80e5),   # Source A
        PressureBoundary(node_id=7, pressure_pa=75e5),   # Source B
    ),
    pressure_outlets  = (
        PressureBoundary(node_id=6, pressure_pa=30e5),   # Gas plant
        PressureBoundary(node_id=8, pressure_pa=40e5),   # Power station
    ),
    components        = pipes,
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

print("=" * 78)
print(f"Case:       {case.name}")
print(f"Converged:  {result.converged}")
print()

NODE_LABELS = {
    1: "Source A (rich)",
    7: "Source B (lean)",
    2: "Ring N2 (A entry)",
    4: "Ring N4 (B entry)",
    3: "Ring N3 (plant tap)",
    5: "Ring N5 (power tap)",
    6: "Gas plant",
    8: "Power station",
}
print(f"{'Node':>5}  {'Label':22}  {'P (bar)':>9}  {'T (°C)':>8}")
print("-" * 50)
for nid in (1, 7, 2, 4, 3, 5, 6, 8):
    P = result.node_pressures_pa[nid]
    T = result.node_temperatures_c.get(nid, float("nan"))
    print(f"{nid:>5}  {NODE_LABELS[nid]:22}  {P/1e5:>9.2f}  {T:>8.2f}")
print()

flows = {cf.label.split(":")[-1]: cf for cf in result.component_flows}

print(f"{'Pipe':11}  {'ṁ (kg/s)':>9}  {'T_in':>7}  {'T_out':>7}  "
      f"{'V gas':>6}  {'L HCliq':>7}  {'W free':>7}  Phase state")
print("-" * 82)
for name, cf in flows.items():
    vf = cf.vapor_fraction  if cf.vapor_fraction  is not None else float("nan")
    lf = cf.liquid_fraction if cf.liquid_fraction is not None else 0.0
    fw = cf.free_water_fraction
    phases = []
    if vf > 1e-3: phases.append("Gas")
    if lf > 1e-3: phases.append("HC liq")
    if fw > 1e-3: phases.append("Free water")
    tin  = f"{cf.temperature_in_c:>7.2f}"  if cf.temperature_in_c  is not None else "    —  "
    tout = f"{cf.temperature_out_c:>7.2f}" if cf.temperature_out_c is not None else "    —  "
    print(f"{name:11}  {cf.mass_flow_kg_per_s:>9.3f}  {tin}  {tout}  "
          f"{vf:>6.3f}  {lf:>7.3f}  {fw:>7.4f}  {' + '.join(phases)}")
print()

# ── How each source is shared between the customers ──────────────────────────
mA, mB = flows["feed_A"].mass_flow_kg_per_s, flows["feed_B"].mass_flow_kg_per_s
print("Source split (mass):")
print(f"  Source A {mA:6.2f} kg/s → plant {flows['ring_north'].mass_flow_kg_per_s / mA:5.1%}, "
      f"power station {flows['ring_west'].mass_flow_kg_per_s / mA:5.1%}")
print(f"  Source B {mB:6.2f} kg/s → plant {flows['ring_east'].mass_flow_kg_per_s / mB:5.1%}, "
      f"power station {flows['ring_south'].mass_flow_kg_per_s / mB:5.1%}")
print()

# ── What each customer receives ──────────────────────────────────────────────
print("Delivered stream (overall mole fractions):")
print(f"  {'Customer':14}  {'ṁ (kg/s)':>9}  " + "  ".join(f"{n[:8]:>8}" for n in COMPONENTS))
for label, key in (("Gas plant", "to_plant"), ("Power station", "to_power")):
    cf = flows[key]
    print(f"  {label:14}  {cf.mass_flow_kg_per_s:>9.2f}  " + "  ".join(f"{v:>8.3f}" for v in cf.zs))
print()

plant, power = flows["to_plant"], flows["to_power"]
print("Delivered phase split (molar):")
print(f"  {'Customer':14}  {'Gas':>6}  {'HC liq':>7}  {'Free water':>10}")
for label, cf in (("Gas plant", plant), ("Power station", power)):
    print(f"  {label:14}  {cf.vapor_fraction:>6.3f}  {cf.liquid_fraction:>7.3f}  {cf.free_water_fraction:>10.4f}")
print()

richer = "gas plant" if plant.zs[0] < power.zs[0] else "power station"
print(f"  • The {richer} sits closer (hydraulically) to Source A and receives the")
print("    richer, wetter blend; the other customer gets a leaner, drier gas.")
print("  • Right-click any pipe in the GUI → Phase Envelope to see where each")
print("    stream sits relative to its bubble and dew curves.")
print()
print(f"Global mass balance error: {result.global_balance.mass_error_pct:.4f} %")
