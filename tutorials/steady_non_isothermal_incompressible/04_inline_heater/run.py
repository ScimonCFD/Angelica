from __future__ import annotations

import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[3] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.cases import build_inline_heater_case
from angelica.core.settings import SolverSettings
from angelica.solvers import SteadyNonIsothermalIncompressibleSolver


def main() -> None:
    case = build_inline_heater_case()
    solver = SteadyNonIsothermalIncompressibleSolver(
        hydraulic_settings=SolverSettings(
            turbulent_iterations=200,
            pressure_relaxation=1.0,
            pressure_correction_abs_tolerance_pa=1e-3,
            pressure_correction_rel_tolerance=1e-6,
        ),
    )
    result = solver.solve(case)

    print(f"Case: {case.name}")
    print(f"Converged:              {result.converged}")
    print(f"Temperature iterations: {len(result.temperature_history)}")
    print()

    # The second component (index 1) is the HeatSource.
    # component_flows carry T_in and T_out at the device boundaries.
    heater = result.component_flows[1]
    mdot = heater.mass_flow_kg_per_s
    T_in  = heater.temperature_in_c
    T_out = heater.temperature_out_c
    cp    = 4182.0
    Q     = 50_000.0

    print(f"Mass flow:  {mdot:.4f} kg/s  ({heater.volumetric_flow_m3_per_h:.3f} m³/h)")
    print(f"T_in  (heater inlet):   {T_in:.4f} °C")
    print(f"T_out (heater outlet):  {T_out:.4f} °C")
    print(f"ΔT (actual):            {T_out - T_in:.4f} K")
    print()

    expected_delta_t = Q / (mdot * cp)
    print(f"Expected ΔT = Q / (ṁ·cp) = {Q:.0f} / ({mdot:.4f} × {cp:.0f}) = {expected_delta_t:.4f} K")
    print(f"Error: {abs((T_out - T_in) - expected_delta_t):.4f} K")


if __name__ == "__main__":
    main()
