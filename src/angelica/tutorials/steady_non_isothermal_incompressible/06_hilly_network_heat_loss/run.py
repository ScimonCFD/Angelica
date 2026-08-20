"""Tutorial 06 — Non-isothermal hilly hot-water network.

Demonstrates that the elevation source term ρ·g·Δz is properly handled by
the non-isothermal solver when pipes ascend and descend by significant amounts.

Network topology
----------------
                         Node 2 (hilltop, +40 m)
                        /          \\
     Node 1 ──trunk──> /            \\──branch A──> Node 3 (outflow, −20 m)
  (600 kPa, 85°C)  ascending +40 m  descending −60 m
                                    \\──branch B──> Node 4 (200 kPa, +10 m)
                                     descending −30 m

Pipe lengths are measured along the pipe axis (actual inclined length):
  trunk:    L=150 m, Δz=+40 m  → inclination ≈ 15.5°
  branch A: L=220 m, Δz=−60 m  → inclination ≈ 15.9°
  branch B: L=120 m, Δz=−30 m  → inclination ≈ 14.5°

Key effects
-----------
  1. Elevation vs friction: ascending the trunk costs ~392 kPa in hydrostatic
     head (ρ·g·40 m ≈ 975 × 9.81 × 40) on top of friction losses.
  2. Gravity-assisted descent: both branches benefit from the downhill head,
     which increases their effective driving pressure.
  3. Outflow BC: node 3 has a prescribed extraction (3.0 kg/s) with no fixed
     pressure — the solver determines its pressure from continuity.
  4. Heat loss: water cools from 85 °C as it travels, and ρ(T) and μ(T) update
     each outer iteration.
"""
from __future__ import annotations

import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[4] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.cases.hilly_hot_water_network import build_hilly_hot_water_network_case
from angelica.solvers import SteadyNonIsothermalIncompressibleSolver

_G = 9.81   # m/s²
_RHO_APPROX = 975.0  # kg/m³ — approximate water density at ~80 °C


def main() -> None:
    case = build_hilly_hot_water_network_case()
    solver = SteadyNonIsothermalIncompressibleSolver()
    result = solver.solve(case)

    print("=" * 62)
    print("Tutorial 06 — Hilly hot-water network (non-isothermal)")
    print("=" * 62)
    print(f"Converged:              {result.converged}")
    print(f"Temperature iterations: {len(result.temperature_history)}")
    print()

    # ── Hydrostatic context ────────────────────────────────────────────────
    h_trunk   = +40.0
    h_branch_a = -60.0
    h_branch_b = -30.0
    print("Expected hydrostatic contributions (ρ·g·Δz):")
    for label, dz in [("trunk (↑)", h_trunk),
                       ("branch A (↓)", h_branch_a),
                       ("branch B (↓)", h_branch_b)]:
        sign = "opposes" if dz > 0 else "assists"
        print(f"  {label:<18}  Δz={dz:+.0f} m  →  "
              f"{abs(_RHO_APPROX*_G*dz)/1e3:.1f} kPa  ({sign} flow)")
    print()

    # ── Nodal results ──────────────────────────────────────────────────────
    print(f"{'Node':<6}  {'P (kPa)':>9}  {'T (°C)':>7}  {'Elevation':>12}")
    elev = {1: 0.0, 2: +40.0, 3: -20.0, 4: +10.0}
    for nid in sorted(result.node_pressures_pa):
        p = result.node_pressures_pa[nid]
        T = result.node_temperatures_c.get(nid, float("nan"))
        z = elev[nid]
        print(f"  {nid:<4}  {p/1e3:>9.2f}  {T:>7.2f}  {z:>+8.1f} m")
    print()

    # ── Component flows ────────────────────────────────────────────────────
    print(f"{'Pipe':<22}  {'ṁ (kg/s)':>10}  {'T_in':>7}  {'T_out':>7}  {'Δz':>8}")
    dz_map = {
        "trunk_ascending":    +40.0,
        "branch_a_descending": -60.0,
        "branch_b_descending": -30.0,
    }
    for cf in result.component_flows:
        cid = cf.label.split(":")[-1]   # strip "Pipe:" prefix
        dz = dz_map.get(cid, 0.0)
        print(f"  {cid:<22}  {cf.mass_flow_kg_per_s:>10.4f}  "
              f"{cf.temperature_in_c:>7.2f}  {cf.temperature_out_c:>7.2f}  "
              f"{dz:>+5.0f} m")
    print()

    # ── Mass balance check ─────────────────────────────────────────────────
    flows = {cf.label.split(":")[-1]: cf.mass_flow_kg_per_s for cf in result.component_flows}
    m_in  = flows.get("trunk_ascending", 0.0)
    m_a   = flows.get("branch_a_descending", 0.0)
    m_b   = flows.get("branch_b_descending", 0.0)
    print(f"Mass balance:")
    print(f"  Trunk in:           {m_in:.4f} kg/s")
    print(f"  Branch A + B out:   {m_a + m_b:.4f} kg/s")
    print(f"  Imbalance:          {abs(m_in - m_a - m_b):.2e} kg/s")
    print()

    # ── Elevation effect interpretation ───────────────────────────────────
    p2 = result.node_pressures_pa[2]
    p3 = result.node_pressures_pa.get(3)
    p4 = result.node_pressures_pa[4]
    print("Elevation effect:")
    print(f"  Node 2 (hilltop):   {p2/1e3:.2f} kPa")
    print(f"  Node 3 (outflow):   {p3/1e3:.2f} kPa  (free — set by continuity)")
    print(f"  Node 4 (outlet):    {p4/1e3:.2f} kPa  (fixed BC)")
    hydro_a = _RHO_APPROX * _G * abs(h_branch_a)
    print(f"  Gravity assist in branch A: {hydro_a/1e3:.1f} kPa")
    print("  Node 3 is BELOW node 1 by 20 m, yet it can sustain pressure")
    print("  because the 60 m descent from node 2 provides significant head.")


if __name__ == "__main__":
    main()
