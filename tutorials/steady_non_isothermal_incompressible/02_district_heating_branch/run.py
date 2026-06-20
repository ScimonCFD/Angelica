from __future__ import annotations

import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[3] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.cases import build_district_heating_branch_case
from angelica.solvers import SteadyNonIsothermalIncompressibleSolver


def main() -> None:
    case = build_district_heating_branch_case()
    solver = SteadyNonIsothermalIncompressibleSolver()
    result = solver.solve(case)

    print(f"Case: {case.name}")
    print(f"Converged:              {result.converged}")
    print(f"Temperature iterations: {len(result.temperature_history)}")
    print()

    print(f"{'Pipe':<25}  {'Mass flow (kg/s)':>16}  {'Vol. flow (m³/h)':>18}")
    labels = ["Main supply (1→2)", "Branch A (2→3)", "Branch B (2→4)"]
    for label, cf in zip(labels, result.component_flows):
        print(
            f"{label:<25}  "
            f"{cf.mass_flow_kg_per_s:>16.4f}  "
            f"{cf.volumetric_flow_m3_per_h:>18.3f}"
        )
    print()

    node_labels = {1: "Supply (source)", 2: "Junction", 3: "Load A (sink)", 4: "Load B (sink)"}
    print(f"{'Node':<20}  {'Pressure (bar)':>14}  {'Temperature (°C)':>16}")
    for nid in sorted(result.node_pressures_pa):
        p_bar = result.node_pressures_pa[nid] / 1e5
        T_c = result.node_temperatures_c[nid]
        label = node_labels.get(nid, str(nid))
        print(f"{label:<20}  {p_bar:>14.4f}  {T_c:>16.2f}")


if __name__ == "__main__":
    main()
