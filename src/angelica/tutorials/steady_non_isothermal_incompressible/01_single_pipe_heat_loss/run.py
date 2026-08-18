from __future__ import annotations

import math
import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[3] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.cases import build_hot_water_pipe_heat_loss_case
from angelica.solvers import SteadyNonIsothermalIncompressibleSolver


def main() -> None:
    case = build_hot_water_pipe_heat_loss_case()
    solver = SteadyNonIsothermalIncompressibleSolver()
    result = solver.solve(case)

    print(f"Case: {case.name}")
    print(f"Converged:               {result.converged}")
    print(f"Temperature iterations:  {len(result.temperature_history)}")
    print()

    mdot = result.component_flows[0].mass_flow_kg_per_s
    vol = result.component_flows[0].volumetric_flow_m3_per_h
    print(f"Mass flow:  {mdot:.4f} kg/s  ({vol:.3f} m³/h)")
    print()

    T_in = result.node_temperatures_c[1]
    T_out = result.node_temperatures_c[2]
    print(f"T_in  (node 1):  {T_in:.2f} °C")
    print(f"T_out (node 2):  {T_out:.2f} °C")
    print(f"ΔT:              {T_in - T_out:.2f} K")
    print()

    # Analytical benchmark (plug-flow with wall heat loss)
    U, D, L, cp, T_amb = 50.0, 0.025, 500.0, 4182.0, 20.0
    T_out_exact = T_amb + (T_in - T_amb) * math.exp(-U * math.pi * D * L / (mdot * cp))
    print(f"Analytical T_out:  {T_out_exact:.2f} °C  (error {abs(T_out - T_out_exact):.2f} K)")


if __name__ == "__main__":
    main()
