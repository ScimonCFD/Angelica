from __future__ import annotations

import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[3] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.cases.crude_oil_pipeline import build_crude_oil_pipeline_case
from angelica.properties import dead_oil_density_kg_per_m3, dead_oil_viscosity_pa_s
from angelica.solvers import SteadyIsothermalIncompressibleSolver

_API = 32.0
_TEMPERATURE_C = 65.0


def main() -> None:
    rho = dead_oil_density_kg_per_m3(_API)
    mu = dead_oil_viscosity_pa_s(_API, _TEMPERATURE_C)

    print(f"Fluid: {_API}°API crude oil at {_TEMPERATURE_C}°C")
    print(f"  Density:   {rho:.2f} kg/m³")
    print(f"  Viscosity: {mu * 1000:.4f} cP  ({mu:.6f} Pa·s)")
    print()

    case = build_crude_oil_pipeline_case()
    solver = SteadyIsothermalIncompressibleSolver()
    result = solver.solve(case)

    print(f"Case: {case.name}")
    print(f"Converged: {result.converged}")
    print(f"Laminar iterations:   {len(result.laminar_metrics)}")
    print(f"Turbulent iterations: {len(result.turbulent_metrics)}")
    print()

    print(f"{'Component':<20}  {'Mass flow (kg/s)':>16}  {'Vol. flow (m³/h)':>18}")
    for comp in result.component_flows:
        print(
            f"{comp.label:<20}  "
            f"{comp.mass_flow_kg_per_s:>16.4f}  "
            f"{comp.volumetric_flow_m3_per_h:>18.4f}"
        )

    print()
    print(f"{'Node':<8}  {'Pressure (Pa)':>14}  {'Pressure (bar)':>14}")
    for nid, p in sorted(result.node_pressures_pa.items()):
        print(f"{nid:<8}  {p:>14.2f}  {p / 1e5:>14.4f}")


if __name__ == "__main__":
    main()
