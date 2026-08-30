"""
Tutorial 02 — Compositional Gas Mixing Junction
================================================
Two gas sources with different compositions feed a common trunk line.
The solver propagates each inlet's mole-fraction vector through the
network and computes a molar-flow-weighted blend at the junction.

Geometry
--------
  Node 1 ──PipeA (5 km, D=0.12 m)──┐
  P=120 bar                          ├── Node 3 ──PipeC (8 km, D=0.18 m)── Node 4
  T=70 °C  rich gas (CH₄ 90%)       │           junction                    P=25 bar
                                     │
  Node 2 ──PipeB (5 km, D=0.10 m)──┘
  P=110 bar
  T=60 °C  leaner gas (CH₄ 60%, C₂H₆ 30%, C₃H₈ 10%)

Fluid
-----
  Three-component system: methane, ethane, propane.

  Source A (rich gas)  : CH₄ 90%,  C₂H₆  8%,  C₃H₈ 2%
  Source B (leaner gas): CH₄ 60%,  C₂H₆ 30%,  C₃H₈ 10%

  At the junction the two streams mix by molar-flow-weighted average.
  The blended composition in PipeC is reported at the end.

What the outer loop does
------------------------
  Each outer iteration:
    1. Propagate inlet compositions → per-pipe PipeState.zs
    2. Flash (P, T, z) at each pipe → density, viscosity, Cp, k
    3. SIMPLE hydraulic solve (inner loop)
    4. Energy equation → update temperatures
    5. Converge when max|Δρ/ρ| and max|ΔT| < tolerance
"""
import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[4] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.core.case import InletCompositionBC, NetworkCase, PressureBoundary, ThermalBoundary
from angelica.core.components import Pipe
from angelica.properties.compositional_fluid import CompositionalFluid
from angelica.solvers import SteadyCompositionalSolver

# ── Fluid definition ──────────────────────────────────────────────────────────
COMPONENTS = ["methane", "ethane", "propane"]

ZS_A = (0.90, 0.08, 0.02)   # Source A — rich gas
ZS_B = (0.60, 0.30, 0.10)   # Source B — leaner gas

# default_zs is used for pipes not yet reached by propagation
fluid = CompositionalFluid(
    components  = COMPONENTS,
    default_zs  = (0.75, 0.19, 0.06),   # rough blend of both sources
)

# ── Network ───────────────────────────────────────────────────────────────────
pipes = [
    Pipe(component_id="pipe_A", start_node=1, end_node=3,
         diameter_m=0.12, length_m=5_000.0, absolute_roughness_m=46e-6,
         heat_transfer_coefficient_w_per_m2k=5.0, ambient_temperature_c=15.0),
    Pipe(component_id="pipe_B", start_node=2, end_node=3,
         diameter_m=0.10, length_m=5_000.0, absolute_roughness_m=46e-6,
         heat_transfer_coefficient_w_per_m2k=5.0, ambient_temperature_c=15.0),
    Pipe(component_id="pipe_C", start_node=3, end_node=4,
         diameter_m=0.18, length_m=8_000.0, absolute_roughness_m=46e-6,
         heat_transfer_coefficient_w_per_m2k=5.0, ambient_temperature_c=15.0),
]

case = NetworkCase(
    name             = "Compositional Gas Mixing Junction",
    fluid_model      = fluid,
    pressure_inlets  = (
        PressureBoundary(node_id=1, pressure_pa=120e5),
        PressureBoundary(node_id=2, pressure_pa=110e5),
    ),
    pressure_outlets = (PressureBoundary(node_id=4, pressure_pa=25e5),),
    components       = tuple(pipes),
    thermal_inlets   = (
        ThermalBoundary(node_id=1, temperature_c=70.0, bc_type="fixed_temperature"),
        ThermalBoundary(node_id=2, temperature_c=60.0, bc_type="fixed_temperature"),
    ),
    inlet_composition_bcs = (
        InletCompositionBC(node_id=1, zs=ZS_A),
        InletCompositionBC(node_id=2, zs=ZS_B),
    ),
)

# ── Solve ─────────────────────────────────────────────────────────────────────
solver = SteadyCompositionalSolver()
result = solver.solve(case)

print("=" * 65)
print(f"Case:              {case.name}")
print(f"Converged:         {result.converged}")
print(f"Outer iterations:  {len(result.density_history)}")
print()

# ── Node results ──────────────────────────────────────────────────────────────
print(f"{'Node':>6}  {'P (bar)':>9}  {'T (°C)':>8}")
for nid in sorted(result.node_pressures_pa):
    P = result.node_pressures_pa[nid]
    T = result.node_temperatures_c[nid]
    print(f"{nid:>6}  {P/1e5:>9.3f}  {T:>8.2f}")
print()

# ── Per-pipe flows ────────────────────────────────────────────────────────────
print(f"{'Pipe':10s}  {'ṁ (kg/s)':>10}  {'Q (m³/h)':>10}")
flows = {}
for cf in result.component_flows:
    label = cf.label.split(":")[-1].strip()
    flows[label] = cf.mass_flow_kg_per_s
    print(f"{label:10s}  {cf.mass_flow_kg_per_s:>10.4f}  {cf.volumetric_flow_m3_per_h:>10.3f}")
print()

# ── Composition at junction: molar-flow-weighted blend ────────────────────────
m_A = flows.get("pipe_A", 0.0)
m_B = flows.get("pipe_B", 0.0)
m_tot = m_A + m_B

if m_tot > 0:
    # Mole fractions are propagated inside the solver's PipeState.zs, but we
    # can also compute the expected blended composition analytically from the
    # individual mass flows and feed compositions.
    # (Assumes similar MWs for the two gases so mass ≈ molar weighting.)
    z_blend = tuple(
        (m_A * z_a + m_B * z_b) / m_tot
        for z_a, z_b in zip(ZS_A, ZS_B)
    )
    print("Blended composition in pipe_C (mass-weighted):")
    for comp, z in zip(COMPONENTS, z_blend):
        print(f"  {comp:10s}: {z*100:.2f} mol%")
    print()

print(f"Global mass balance error: {result.global_balance.mass_error_pct:.4f} %")
