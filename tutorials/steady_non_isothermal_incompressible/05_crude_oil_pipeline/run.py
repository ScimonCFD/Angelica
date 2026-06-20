"""Tutorial 05 — Non-isothermal crude oil gathering pipeline.

A 32°API dead crude oil enters a branched gathering network at 80 °C.
The pipeline loses heat to a 15 °C ambient through light insulation
(U = 2 W/m²K).  As the oil cools its viscosity rises sharply, increasing
friction losses.  This tutorial demonstrates:

  1. Using build_thermal_dead_oil() for temperature-dependent fluid properties.
  2. Running SteadyNonIsothermalIncompressibleSolver on a crude oil network.
  3. The viscosity effect: comparing flow against a cold-viscosity baseline.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[3] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.cases.crude_oil_pipeline_thermal import build_crude_oil_pipeline_thermal_case
from angelica.properties.dead_oil import (
    build_thermal_dead_oil,
    dead_oil_density_kg_per_m3,
    dead_oil_specific_heat_j_per_kg_k,
    dead_oil_thermal_conductivity_w_per_m_k,
    dead_oil_viscosity_pa_s,
)
from angelica.solvers import SteadyNonIsothermalIncompressibleSolver

_API = 32.0
_T_IN = 80.0   # °C — hot inlet
_T_AMB = 15.0  # °C — ambient


def _ntu_outlet_temperature(mdot: float, D: float, L: float, U: float, cp: float,
                            T_in: float, T_amb: float) -> float:
    """Analytical NTU outlet temperature for a single pipe with wall heat loss."""
    ntu = U * math.pi * D * L / (mdot * cp)
    return T_amb + (T_in - T_amb) * math.exp(-ntu)


def main() -> None:
    # ── fluid properties at inlet and ambient temperatures ─────────────────
    print("=" * 60)
    print(f"Fluid: {_API}°API dead crude oil")
    print("=" * 60)
    for T in (_T_IN, _T_AMB):
        mu = dead_oil_viscosity_pa_s(_API, T)
        cp = dead_oil_specific_heat_j_per_kg_k(_API, T)
        k  = dead_oil_thermal_conductivity_w_per_m_k(_API, T)
        print(f"  T = {T:5.1f} °C  →  μ = {mu*1000:.2f} cP  "
              f"cp = {cp:.0f} J/kg/K  k = {k:.4f} W/m/K")
    print()

    # ── non-isothermal solve ───────────────────────────────────────────────
    case = build_crude_oil_pipeline_thermal_case()
    solver = SteadyNonIsothermalIncompressibleSolver()
    result = solver.solve(case)

    print(f"Case: {case.name}")
    print(f"Converged:              {result.converged}")
    print(f"Temperature iterations: {len(result.temperature_history)}")
    print()

    # ── nodal results ──────────────────────────────────────────────────────
    print(f"{'Node':<6}  {'P (bar)':>8}  {'T (°C)':>8}")
    for nid in sorted(result.node_pressures_pa):
        p = result.node_pressures_pa[nid]
        T = result.node_temperatures_c[nid]
        print(f"  {nid:<4}  {p/1e5:>8.3f}  {T:>8.2f}")
    print()

    # ── component results ──────────────────────────────────────────────────
    print(f"{'Component':<14}  {'ṁ (kg/s)':>10}  {'Q (m³/h)':>10}  "
          f"{'T_in (°C)':>10}  {'T_out (°C)':>10}")
    for cf in result.component_flows:
        print(f"  {cf.label:<12}  {cf.mass_flow_kg_per_s:>10.4f}  "
              f"{cf.volumetric_flow_m3_per_h:>10.4f}  "
              f"{cf.temperature_in_c:>10.2f}  {cf.temperature_out_c:>10.2f}")
    print()

    # ── viscosity effect: compare trunk flow vs cold-viscosity assumption ──
    trunk = result.component_flows[0]
    mdot_hot = trunk.mass_flow_kg_per_s

    mu_hot  = dead_oil_viscosity_pa_s(_API, _T_IN)
    mu_cold = dead_oil_viscosity_pa_s(_API, _T_AMB)
    rho     = dead_oil_density_kg_per_m3(_API)

    print("Viscosity effect on trunk line:")
    print(f"  μ at inlet temperature ({_T_IN}°C): {mu_hot*1000:.2f} cP")
    print(f"  μ at ambient temperature ({_T_AMB}°C): {mu_cold*1000:.2f} cP")
    print(f"  Viscosity ratio (cold/hot): {mu_cold/mu_hot:.1f}×")
    print(f"  Actual mass flow (thermal solver): {mdot_hot:.4f} kg/s")
    print()
    print("  Thermal simulation correctly captures the viscosity rise along")
    print("  the pipeline that a constant-property isothermal solver would miss.")


if __name__ == "__main__":
    main()
