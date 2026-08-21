"""
Tutorial 03 — Compositional Looped Gas Network
===============================================
Two gas sources with different compositions feed a looped gathering
network.  The loop creates a hydraulic cycle — the SIMPLE solver
distributes flow to satisfy mass balance at every node and Kirchhoff's
pressure law around the loop.

Geometry
--------

  Node 1 ──PipeA (5 km, D=0.12 m)──→ Node 3 ──PipeC (8 km, D=0.15 m)──→ Node 5
  P=120 bar,                             │                                   (outlet,
  T=70 °C                               PipeD (3 km, D=0.10 m)              P=20 bar)
  Rich gas:                              │ ← loop pipe                         ↑
    CH₄  90%                           Node 4 ──PipeE (8 km, D=0.15 m)────┘
    C₂H₆  8%                             ↑
    C₃H₈  2%                           PipeB (5 km, D=0.12 m)
                                          │
  Node 2  P=110 bar, T=60 °C
  Lean gas:  CH₄ 60%, C₂H₆ 30%, C₃H₈ 10%

Key things to observe
---------------------
  1. PipeD direction — the loop pipe carries gas from the
     higher-pressure junction to the lower-pressure one.  Which way
     it flows is not prescribed; it emerges from the hydraulic solution.

  2. Composition at Node 4 — it receives both the lean gas from PipeB
     and whatever PipeD brings from Node 3.  The blend depends on the
     mass-flow ratio, which the solver computes.

  3. Delivered composition at Node 5 — a mixture of both sources,
     arriving via PipeC and PipeE with different weights.
"""
import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[4] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.cases.looped_gas_gathering import build_looped_gas_gathering_case
from angelica.solvers import SteadyCompositionalSolver

# ── Solve ─────────────────────────────────────────────────────────────────────
case   = build_looped_gas_gathering_case()
solver = SteadyCompositionalSolver()
result = solver.solve(case)

print("=" * 65)
print(f"Case:              {case.name}")
print(f"Converged:         {result.converged}")
print(f"Outer iterations:  {len(result.density_history)}")
print()

# ── Node pressures and temperatures ──────────────────────────────────────────
labels = {1: "Source A",  2: "Source B",
          3: "Junction 3", 4: "Junction 4", 5: "Outlet"}
print(f"{'Node':>6}  {'Label':14s}  {'P (bar)':>9}  {'T (°C)':>8}")
for nid in sorted(result.node_pressures_pa):
    P = result.node_pressures_pa[nid]
    T = result.node_temperatures_c[nid]
    print(f"{nid:>6}  {labels[nid]:14s}  {P/1e5:>9.3f}  {T:>8.2f}")
print()

# ── Per-pipe flows ────────────────────────────────────────────────────────────
flows = {}
for cf in result.component_flows:
    label = cf.label.split(":")[-1].strip()
    flows[label] = cf.mass_flow_kg_per_s

loop_pipe = "pipe_D"
m_D = flows.get(loop_pipe, 0.0)

print(f"{'Pipe':10s}  {'ṁ (kg/s)':>10}  {'Q (m³/h)':>10}  Direction")
for cf in result.component_flows:
    label = cf.label.split(":")[-1].strip()
    m     = cf.mass_flow_kg_per_s
    q     = cf.volumetric_flow_m3_per_h
    if label == loop_pipe:
        direction = "3 → 4" if m >= 0 else "4 → 3  ← reversed!"
        direction = f"** {direction} **"
    else:
        direction = ""
    print(f"{label:10s}  {m:>10.4f}  {q:>10.3f}  {direction}")
print()

# ── Junction compositions (analytical from flow field) ────────────────────────
# Inlet compositions (mole fractions)
ZS_A = (0.90, 0.08, 0.02)   # Source A
ZS_B = (0.60, 0.30, 0.10)   # Source B
COMPS = ["methane", "ethane", "propane"]

m_A = flows.get("pipe_A", 0.0)
m_B = flows.get("pipe_B", 0.0)

def blend(streams):
    """Mass-weighted mole-fraction blend of [(mdot, zs), ...]."""
    total = sum(m for m, _ in streams if m > 0)
    if total <= 0:
        return None
    return tuple(
        sum(m * z[i] for m, z in streams if m > 0) / total
        for i in range(len(streams[0][1]))
    )

# Node 3: fed by pipe_A; pipe_D may bring additional gas if it flows 4→3
if m_D >= 0:
    z3 = ZS_A                              # only PipeA feeds Node 3
    z4 = blend([(m_B, ZS_B), (m_D, ZS_A)])  # PipeB + PipeD from Node 3
else:
    # PipeD is reversed: Node 4 feeds Node 3
    z4 = ZS_B                              # only PipeB feeds Node 4
    z3 = blend([(m_A, ZS_A), (-m_D, ZS_B)])  # PipeA + PipeD from Node 4

# Node 5: fed by pipe_C (from Node 3) and pipe_E (from Node 4)
m_C = flows.get("pipe_C", 0.0)
m_E = flows.get("pipe_E", 0.0)
z5  = blend([(m_C, z3), (m_E, z4)])

print("Compositions at junction and outlet nodes (mol%):")
header = f"{'Component':12s}  {'Node 3':>10}  {'Node 4':>10}  {'Node 5 (out)':>12}"
print(header)
print("-" * len(header))
for i, comp in enumerate(COMPS):
    v3 = z3[i] * 100 if z3 else 0
    v4 = z4[i] * 100 if z4 else 0
    v5 = z5[i] * 100 if z5 else 0
    print(f"{comp:12s}  {v3:>10.2f}  {v4:>10.2f}  {v5:>12.2f}")
print()

print(f"Global mass balance error: {result.global_balance.mass_error_pct:.4f} %")
