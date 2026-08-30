from __future__ import annotations

import sys
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[3] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.cases import build_looped_network_heat_loss_case
from angelica.core.settings import SolverSettings
from angelica.solvers import NonIsothermalSolverSettings, SteadyNonIsothermalIncompressibleSolver


def main() -> None:
    case = build_looped_network_heat_loss_case()

    # temperature_relaxation=0.5 slows the outer loop so convergence is visible
    # on the plot (~12 iterations instead of 3–4).
    solver = SteadyNonIsothermalIncompressibleSolver(
        hydraulic_settings=SolverSettings(
            turbulent_iterations=200,
            pressure_relaxation=0.7,
            pressure_correction_abs_tolerance_pa=1e-3,
            pressure_correction_rel_tolerance=1e-8,
        ),
        non_isothermal_settings=NonIsothermalSolverSettings(
            temperature_relaxation=0.5,
            max_temperature_iterations=50,
        ),
    )
    result = solver.solve(case)

    print(f"Case: {case.name}")
    print(f"Converged:              {result.converged}")
    print(f"Temperature iterations: {len(result.temperature_history)}")
    print(f"Hyd. iterations:        {len(result.turbulent_metrics)}")
    print()

    node_labels = {
        1: "Source",
        2: "Junction A",
        3: "Junction B",
        4: "Junction C (mixing)",
        5: "Sink",
    }
    print(f"{'Node':<28}  {'P (bar)':>8}  {'T (°C)':>8}")
    for nid in sorted(result.node_pressures_pa):
        p = result.node_pressures_pa[nid] / 1e5
        T = result.node_temperatures_c[nid]
        print(f"  {node_labels[nid]:<26}  {p:>8.3f}  {T:>8.2f}")
    print()

    pipe_labels = [
        "Feed  1→2 (D50, 100m)",
        "Upper 2→3 (D40, 400m)",
        "Upper 3→4 (D40, 300m)",
        "Bypass 2→4 (D35, 400m)",
        "Exit  4→5 (D50, 100m)",
    ]
    print(f"{'Pipe':<28}  {'kg/s':>8}  {'m³/h':>8}")
    for label, cf in zip(pipe_labels, result.component_flows):
        print(f"  {label:<26}  {cf.mass_flow_kg_per_s:>8.4f}  {cf.volumetric_flow_m3_per_h:>8.3f}")
    print()

    print("Temperature convergence history (max |ΔT| per outer iteration):")
    for i, delta in enumerate(result.temperature_history, 1):
        bar = "█" * max(1, int(40 * delta / result.temperature_history[0]))
        print(f"  iter {i:2d}: {delta:8.4f} K  {bar}")


if __name__ == "__main__":
    main()
