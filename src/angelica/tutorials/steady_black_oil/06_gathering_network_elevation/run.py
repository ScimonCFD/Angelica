"""Tutorial 06 — Black-oil gathering network with topographic elevation.

Demonstrates the elevation source term ρ_mix·g·Δz in a three-phase system
where the two wells are at very different elevations.  Well A flows downhill
(gravity assists), Well B flows uphill (gravity opposes) to reach the manifold.

Network topology
----------------
  Node 1 (Well A, +150 m)                Node 2 (Well B, −100 m)
  P=6 MPa, T=65°C, 32°API                P=9 MPa, T=60°C, 28°API
       \\                                      /
        \\  descending −150 m    ascending +100 m
         \\                                  /
          ──────────> Node 3 (Manifold, 0 m) <──────────
                           |
                      descending −30 m
                           |
                      Node 4 (Separator, −30 m, P=3 MPa)

Elevation contributions (ρ_mix ≈ 700 kg/m³):
  Well A flowline  (Δz=−150 m): ρ·g·150 ≈ 1.03 MPa  assists flow
  Well B flowline  (Δz=+100 m): ρ·g·100 ≈ 0.69 MPa  opposes flow

Despite Well A having a lower wellhead pressure (6 MPa vs 9 MPa), it can
still deliver significant flow because it gains ~1 MPa from gravity.
Well B needs higher wellhead pressure to overcome its uphill journey.
"""
from __future__ import annotations

import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[4] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.cases.black_oil_gathering_elevation import build_black_oil_gathering_elevation_case
from angelica.properties.black_oil import bubble_point_pa
from angelica.solvers import SteadyBlackOilSolver

_G = 9.81   # m/s²

# ── Fluid parameters for display ──────────────────────────────────────────────
_API_A, _GAS_A, _GOR_A, _T_A = 32.0, 0.65, 25.0, 65.0
_API_B, _GAS_B, _GOR_B, _T_B = 28.0, 0.68, 15.0, 60.0


def main() -> None:
    Pb_A = bubble_point_pa(_GOR_A, _GAS_A, _API_A, _T_A)
    Pb_B = bubble_point_pa(_GOR_B, _GAS_B, _API_B, _T_B)

    print("=" * 64)
    print("Tutorial 06 — Black-oil gathering with elevation (two wells)")
    print("=" * 64)
    print(f"Well A ({_API_A}°API, GOR={_GOR_A}): bubble point = {Pb_A/1e6:.2f} MPa at {_T_A}°C")
    print(f"Well B ({_API_B}°API, GOR={_GOR_B}): bubble point = {Pb_B/1e6:.2f} MPa at {_T_B}°C")
    print()

    # ── Elevation gravity estimate ─────────────────────────────────────────
    rho_approx = 700.0  # kg/m³ — typical live-oil mixture at ~5 MPa
    for label, dz in [("Well A flowline (↓ 150 m)", -150.0),
                       ("Well B flowline (↑ 100 m)", +100.0)]:
        effect = rho_approx * _G * abs(dz) / 1e6
        sign = "assists" if dz < 0 else "opposes"
        print(f"  ρ_mix·g·|Δz| for {label:<28}: {effect:.2f} MPa  ({sign} flow)")
    print()

    # ── Solve ─────────────────────────────────────────────────────────────
    case = build_black_oil_gathering_elevation_case()
    solver = SteadyBlackOilSolver()
    result = solver.solve(case)

    print(f"Converged:          {result.converged}")
    print(f"Density iterations: {len(result.density_history)}")
    print()

    # ── Nodal results ─────────────────────────────────────────────────────
    elev = {1: +150.0, 2: -100.0, 3: 0.0, 4: -30.0}
    print(f"{'Node':<6}  {'P (MPa)':>8}  {'T (°C)':>7}  {'Elevation':>12}")
    for nid in sorted(result.node_pressures_pa):
        p = result.node_pressures_pa[nid]
        T = result.node_temperatures_c.get(nid, float("nan"))
        z = elev[nid]
        print(f"  {nid:<4}  {p/1e6:>8.4f}  {T:>7.2f}  {z:>+8.1f} m")
    print()

    # ── Component flows ───────────────────────────────────────────────────
    print(f"{'Pipe':<24}  {'ṁ (kg/s)':>10}  {'T_in':>7}  {'T_out':>7}  {'Δz':>8}")
    dz_map = {
        "well_a_flowline":        -150.0,
        "well_b_flowline":        +100.0,
        "manifold_to_separator":   -30.0,
    }
    for cf in result.component_flows:
        cid = cf.label.split(":")[-1]
        dz = dz_map.get(cid, 0.0)
        print(f"  {cid:<24}  {cf.mass_flow_kg_per_s:>10.4f}  "
              f"{cf.temperature_in_c:>7.2f}  {cf.temperature_out_c:>7.2f}  "
              f"{dz:>+5.0f} m")
    print()

    # ── Mass balance ──────────────────────────────────────────────────────
    flows = {cf.label.split(":")[-1]: cf.mass_flow_kg_per_s for cf in result.component_flows}
    m_a   = flows.get("well_a_flowline", 0.0)
    m_b   = flows.get("well_b_flowline", 0.0)
    m_sep = flows.get("manifold_to_separator", 0.0)
    print(f"Mass balance:")
    print(f"  Well A + Well B in:  {m_a + m_b:.4f} kg/s")
    print(f"  Separator out:       {m_sep:.4f} kg/s")
    print(f"  Imbalance:           {abs(m_a + m_b - m_sep):.2e} kg/s")
    print()

    # ── Flow split interpretation ─────────────────────────────────────────
    total = m_a + m_b
    if total > 0:
        frac_a = m_a / total * 100.0
        frac_b = m_b / total * 100.0
        print("Flow split at manifold:")
        print(f"  Well A (6 MPa, ↓ 150 m): {m_a:.4f} kg/s  ({frac_a:.1f}%)")
        print(f"  Well B (9 MPa, ↑ 100 m): {m_b:.4f} kg/s  ({frac_b:.1f}%)")
        print("  Despite lower wellhead pressure, Well A compensates via gravity.")


if __name__ == "__main__":
    main()
