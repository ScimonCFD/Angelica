from __future__ import annotations

import math
import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[3] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.cases.three_reservoir_junction import build_three_reservoir_junction_case
from angelica.closures import ColebrookPipeCorrelation
from angelica.core.settings import SolverSettings
from angelica.solvers import SteadyIsothermalIncompressibleSolver

# ── Fluid and geometry constants ────────────────────────────────────────────
_RHO = 998.25   # kg/m³
_MU  = 0.001    # Pa·s
_G   = 9.81     # m/s²

# Pipe definitions: (diameter_m, length_m, roughness_m)
_PIPES = {
    "A→J": (0.200, 1000.0, 0.0001),
    "J→B": (0.150,  800.0, 0.0001),
    "J→C": (0.100,  600.0, 0.0001),
}

# Reservoir heads (m)
_H_A = 30.0
_H_B = 20.0
_H_C =  5.0


# ── Analytical helpers ──────────────────────────────────────────────────────

def _colebrook_white(Re: float, D: float, k: float, tol: float = 1e-10) -> float:
    """Solve Colebrook-White for Darcy friction factor f.

    1/√f = -2 log10(ε/(3.7 D) + 2.51/(Re √f))
    """
    if Re < 1e-12:
        return 0.0
    # Swamee-Jain initial guess
    f = 0.25 / (math.log10(k / (3.7 * D) + 5.74 / Re ** 0.9)) ** 2
    for _ in range(100):
        rhs = -2.0 * math.log10(k / (3.7 * D) + 2.51 / (Re * math.sqrt(f)))
        f_new = 1.0 / rhs ** 2
        if abs(f_new - f) < tol:
            return f_new
        f = f_new
    return f


def _flow_for_head_diff(dH: float, D: float, L: float, k: float) -> float:
    """Return signed volumetric flow (m³/s) for head difference dH (upstream minus downstream).

    Uses Darcy-Weisbach + Colebrook-White iteratively.
    dH > 0  →  flow is positive (in direction of head drop).
    """
    if abs(dH) < 1e-14:
        return 0.0
    dp = abs(dH) * _RHO * _G   # Pa
    A  = math.pi * D ** 2 / 4.0

    # Initial guess via Swamee-Jain
    f = 0.02
    for _ in range(200):
        v = math.sqrt(2.0 * dp * D / (_RHO * f * L))
        Q = v * A
        Re = _RHO * v * D / _MU
        if Re < 1e-12:
            return 0.0
        f_new = _colebrook_white(Re, D, k)
        if abs(f_new - f) < 1e-12:
            f = f_new
            break
        f = f_new

    v = math.sqrt(2.0 * dp * D / (_RHO * f * L))
    Q = v * A
    return math.copysign(Q, dH)


def _mass_balance_at_J(H_J: float) -> float:
    """Net mass flow into junction J for a given junction head H_J.

    Positive = more flow in than out → H_J should decrease.
    Returns Q_AJ - Q_JB - Q_JC  (all in m³/s, positive means toward J).
    """
    D_AJ, L_AJ, k_AJ = _PIPES["A→J"]
    D_JB, L_JB, k_JB = _PIPES["J→B"]
    D_JC, L_JC, k_JC = _PIPES["J→C"]

    Q_AJ = _flow_for_head_diff(_H_A - H_J, D_AJ, L_AJ, k_AJ)
    Q_JB = _flow_for_head_diff(H_J - _H_B, D_JB, L_JB, k_JB)
    Q_JC = _flow_for_head_diff(H_J - _H_C, D_JC, L_JC, k_JC)
    return Q_AJ - Q_JB - Q_JC


def solve_analytical() -> dict:
    """Find H_J by bisection, then compute pipe flows."""
    # Bisect on (_H_B, _H_A) — J must lie between B and A heads
    lo, hi = _H_B + 1e-9, _H_A - 1e-9
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if _mass_balance_at_J(mid) > 0:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-12:
            break
    H_J = 0.5 * (lo + hi)

    D_AJ, L_AJ, k_AJ = _PIPES["A→J"]
    D_JB, L_JB, k_JB = _PIPES["J→B"]
    D_JC, L_JC, k_JC = _PIPES["J→C"]

    Q_AJ = _flow_for_head_diff(_H_A - H_J, D_AJ, L_AJ, k_AJ)
    Q_JB = _flow_for_head_diff(H_J - _H_B, D_JB, L_JB, k_JB)
    Q_JC = _flow_for_head_diff(H_J - _H_C, D_JC, L_JC, k_JC)

    return {
        "H_J":  H_J,
        "Q_AJ": Q_AJ * 3600.0,   # → m³/h
        "Q_JB": Q_JB * 3600.0,
        "Q_JC": Q_JC * 3600.0,
    }


# ── Angelica solve ───────────────────────────────────────────────────────────

def solve_angelica() -> tuple:
    case = build_three_reservoir_junction_case()
    solver = SteadyIsothermalIncompressibleSolver(
        settings=SolverSettings(
            turbulent_iterations=300,
            pressure_relaxation=0.5,
            pressure_correction_abs_tolerance_pa=1e-3,
            pressure_correction_rel_tolerance=1e-8,
        ),
        turbulent_pipe_correlation=ColebrookPipeCorrelation(),
    )
    result = solver.solve(case)

    # Node 2 is the junction J
    H_J = result.node_pressures_pa[2] / (_RHO * _G)
    flows = result.component_flows   # A→J, J→B, J→C
    Q_AJ = flows[0].volumetric_flow_m3_per_h
    Q_JB = flows[1].volumetric_flow_m3_per_h
    Q_JC = flows[2].volumetric_flow_m3_per_h

    return result, {
        "H_J":  H_J,
        "Q_AJ": Q_AJ,
        "Q_JB": Q_JB,
        "Q_JC": Q_JC,
    }


# ── Main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    ref = solve_analytical()
    result, ang = solve_angelica()

    print(f"Converged: {result.converged}")
    print(f"Laminar iterations:   {len(result.laminar_metrics)}")
    print(f"Turbulent iterations: {len(result.turbulent_metrics)}")
    print()

    col = (26, 12, 12, 10)
    hdr = f"{'Quantity':<{col[0]}}  {'Analytical':>{col[1]}}  {'Angelica':>{col[2]}}  {'Delta':>{col[3]}}"
    sep = "-" * (sum(col) + 3 * 2)
    print(hdr)
    print(sep)

    rows = [
        ("H_J (m)",          ref["H_J"],  ang["H_J"],  ""),
        ("Flow A→J (m³/h)",  ref["Q_AJ"], ang["Q_AJ"], ""),
        ("Flow J→B (m³/h)",  ref["Q_JB"], ang["Q_JB"], ""),
        ("Flow J→C (m³/h)",  ref["Q_JC"], ang["Q_JC"], ""),
    ]

    for label, a_val, b_val, _ in rows:
        delta = b_val - a_val
        print(
            f"{label:<{col[0]}}  "
            f"{a_val:>{col[1]}.5f}  "
            f"{b_val:>{col[2]}.5f}  "
            f"{delta:>{col[3]}.5f}"
        )


if __name__ == "__main__":
    main()
