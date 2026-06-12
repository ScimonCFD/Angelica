from __future__ import annotations

import math
import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[3] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.cases.laminar_parallel_pipes import build_laminar_parallel_pipes_case
from angelica.solvers import SteadyIsothermalIncompressibleSolver

# ── Fluid and geometry constants ────────────────────────────────────────────
_RHO   = 880.0   # kg/m³
_MU    = 0.1     # Pa·s
_DP    = 50000.0 # Pa  (P_in - P_out)

# Pipe definitions: (label, D m, L m)
_PIPES = [
    ("Pipe 1", 0.020, 10.0),
    ("Pipe 2", 0.015,  8.0),
    ("Pipe 3", 0.025, 12.0),
]


# ── Analytical Poiseuille solution ──────────────────────────────────────────

def poiseuille_flow_m3_per_s(D: float, L: float, dp: float, mu: float) -> float:
    """Hagen-Poiseuille volumetric flow rate: Q = π D⁴ ΔP / (128 μ L)."""
    return math.pi * D**4 * dp / (128.0 * mu * L)


def reynolds(D: float, Q_m3s: float) -> float:
    """Re = ρ v D / μ  =  ρ Q D / (A μ)  =  4 ρ Q / (π D μ)."""
    return 4.0 * _RHO * Q_m3s / (math.pi * D * _MU)


# ── Angelica solve ───────────────────────────────────────────────────────────

def solve_angelica():
    case = build_laminar_parallel_pipes_case()
    result = SteadyIsothermalIncompressibleSolver().solve(case)
    return result


# ── Main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    # Analytical reference
    analytical = []
    for label, D, L in _PIPES:
        Q_m3s = poiseuille_flow_m3_per_s(D, L, _DP, _MU)
        Q_m3h = Q_m3s * 3600.0
        Re = reynolds(D, Q_m3s)
        analytical.append((label, Q_m3h, Re))

    # Angelica
    result = solve_angelica()
    ang_flows = [comp.volumetric_flow_m3_per_h for comp in result.component_flows]

    print(f"Converged: {result.converged}")
    print(f"Laminar iterations:   {len(result.laminar_metrics)}")
    print(f"Turbulent iterations: {len(result.turbulent_metrics)}")
    print()

    # Verify laminar regime
    all_laminar = all(Re < 2300.0 for _, _, Re in analytical)
    max_re = max(Re for _, _, Re in analytical)
    print(f"All Re < 2300 (laminar regime confirmed): {all_laminar}  (max Re = {max_re:.1f})")
    print()

    # Comparison table
    col = (8, 20, 20, 16, 10)
    hdr = (
        f"{'Pipe':<{col[0]}}  "
        f"{'Analytical (m³/h)':>{col[1]}}  "
        f"{'Angelica (m³/h)':>{col[2]}}  "
        f"{'Delta (m³/h)':>{col[3]}}  "
        f"{'Re':>{col[4]}}"
    )
    sep = "-" * (sum(col) + 4 * 2)
    print(hdr)
    print(sep)

    max_err_pct = 0.0
    for i, ((label, Q_ref, Re_val), Q_ang) in enumerate(zip(analytical, ang_flows)):
        delta = Q_ang - Q_ref
        err_pct = abs(delta) / Q_ref * 100.0
        max_err_pct = max(max_err_pct, err_pct)
        print(
            f"{label:<{col[0]}}  "
            f"{Q_ref:>{col[1]}.6f}  "
            f"{Q_ang:>{col[2]}.6f}  "
            f"{delta:>{col[3]}.6f}  "
            f"{Re_val:>{col[4]}.1f}"
        )

    print()
    print(f"Max relative error: {max_err_pct:.4f} %")


if __name__ == "__main__":
    main()
