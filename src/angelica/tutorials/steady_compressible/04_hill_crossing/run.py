"""Tutorial 04 — Compressible gas pipeline hill crossing.

Demonstrates that the elevation source term ρ·g·Δz is correctly handled
in the compressible solver where gas density varies with pressure.

Network topology
----------------
                    Node 2 (hilltop, +500 m)
                   /           \\
  Node 1 ──────> /               \\──branch A──> Node 3 (+200 m, 5.5 MPa outlet)
 (8 MPa, 15°C)  ascending       descending −300 m
                +500 m            \\──branch B──> Node 4 (−50 m, outflow 8 kg/s)
                                   descending −550 m

Pipe lengths are measured along the pipe axis:
  trunk:    L=600 m, Δz=+500 m  → inclination ≈ 56°
  branch A: L=380 m, Δz=−300 m  → inclination ≈ 52°
  branch B: L=640 m, Δz=−550 m  → inclination ≈ 59°

(These are steep mountain crossings — representative of challenging terrain.)

Key effects
-----------
  1. Gas elevation correction: for methane at 8 MPa, 15 °C, ρ ≈ 53 kg/m³.
     ρ·g·500 m ≈ 260 kPa — a measurable fraction (~3%) of the inlet pressure.
  2. Density-updated gravity: as P drops and gas expands, ρ decreases, so the
     hydrostatic correction is recalculated each outer iteration.
  3. Branch B descends past the supply station (Δz = −550 m total), so the
     descent partially recovers the head lost in the ascent.
  4. Outflow BC: node 4 has prescribed ṁ_out = 8 kg/s with free pressure.
"""
from __future__ import annotations

import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[4] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.cases.gas_pipeline_hill_crossing import build_gas_pipeline_hill_crossing_case
from angelica.properties.eos import IdealGasEOS
from angelica.solvers import SteadyCompressibleSolver

_M_METHANE = 0.016043   # kg/mol
_G = 9.81               # m/s²


def rho_ideal(pressure_pa: float, temperature_c: float = 15.0) -> float:
    eos = IdealGasEOS(molecular_weight_kg_per_mol=_M_METHANE)
    return eos.density(pressure_pa, temperature_c)


def main() -> None:
    case = build_gas_pipeline_hill_crossing_case()
    solver = SteadyCompressibleSolver()
    result = solver.solve(case)

    print("=" * 62)
    print("Tutorial 04 — Gas pipeline hill crossing (compressible)")
    print("=" * 62)
    print(f"Converged:          {result.converged}")
    print(f"Density iterations: {len(result.density_history)}")
    print()

    # ── Elevation-density context ──────────────────────────────────────────
    p_in = result.node_pressures_pa[1]
    rho_in = rho_ideal(p_in)
    print(f"Inlet density at {p_in/1e6:.1f} MPa, 15°C: {rho_in:.2f} kg/m³")
    for label, dz in [("trunk (↑ 500 m)", +500.0),
                       ("branch A (↓ 300 m)", -300.0),
                       ("branch B (↓ 550 m)", -550.0)]:
        corr_kpa = rho_in * _G * abs(dz) / 1e3
        sign = "opposes" if dz > 0 else "assists"
        print(f"  {label:<22}  ρ·g·|Δz| ≈ {corr_kpa:.1f} kPa  ({sign} flow)")
    print()

    # ── Nodal pressures and densities ────────────────────────────────────
    elev = {1: 0, 2: +500, 3: +200, 4: -50}
    print(f"{'Node':<6}  {'P (MPa)':>8}  {'ρ (kg/m³)':>10}  {'Elevation':>12}")
    for nid in sorted(result.node_pressures_pa):
        p = result.node_pressures_pa[nid]
        rho = rho_ideal(p)
        z = elev[nid]
        print(f"  {nid:<4}  {p/1e6:>8.4f}  {rho:>10.3f}  {z:>+8.0f} m")
    print()

    # ── Component flows ──────────────────────────────────────────────────
    print(f"{'Pipe':<26}  {'ṁ (kg/s)':>10}  {'Q (m³/h)':>10}  {'Δz':>8}")
    dz_map = {
        "trunk_ascending":          +500.0,
        "branch_a_partial_descent": -300.0,
        "branch_b_full_descent":    -550.0,
    }
    for cf in result.component_flows:
        cid = cf.label.split(":")[-1]
        dz = dz_map.get(cid, 0.0)
        print(f"  {cid:<26}  {cf.mass_flow_kg_per_s:>10.4f}  "
              f"{cf.volumetric_flow_m3_per_h:>10.3f}  {dz:>+5.0f} m")
    print()

    # ── Mass balance ─────────────────────────────────────────────────────
    flows = {cf.label.split(":")[-1]: cf.mass_flow_kg_per_s for cf in result.component_flows}
    m_trunk = flows.get("trunk_ascending", 0.0)
    m_a     = flows.get("branch_a_partial_descent", 0.0)
    m_b     = flows.get("branch_b_full_descent", 0.0)
    print(f"Mass balance:")
    print(f"  Trunk in:            {m_trunk:.4f} kg/s")
    print(f"  Branch A + B out:    {m_a + m_b:.4f} kg/s")
    print(f"  Imbalance:           {abs(m_trunk - m_a - m_b):.2e} kg/s")
    print()

    # ── Hill-crossing note ───────────────────────────────────────────────
    p2 = result.node_pressures_pa[2]
    p4 = result.node_pressures_pa.get(4)
    print("Hill-crossing observation:")
    print(f"  Pressure at hilltop (node 2): {p2/1e6:.4f} MPa")
    print(f"  Pressure at valley  (node 4): {p4/1e6:.4f} MPa  (free — set by continuity)")
    print("  Branch B descends 550 m past the inlet level, so the gravity")
    print("  term in that pipe partially recovers the head lost climbing the hill.")


if __name__ == "__main__":
    main()
