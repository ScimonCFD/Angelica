from __future__ import annotations

import sys
import unittest
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[1] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.cases import (
    build_steady_water_network_aggressive_elevation_case,
    build_steady_water_network_aggressive_elevation_outlet_flow_case,
    build_steady_water_network_case,
    build_steady_water_network_inlet_flow_boundary_case,
    build_steady_water_network_no_fittings_case,
    build_steady_water_network_two_flow_boundaries_case,
    build_three_reservoir_junction_case,
)
from angelica.closures import ColebrookPipeCorrelation
from angelica.core.settings import SolverSettings
from angelica.solvers import SteadyIsothermalIncompressibleSolver


def solve(case, settings: SolverSettings):
    solver = SteadyIsothermalIncompressibleSolver(settings)
    return solver.solve(case)


class TutorialSuiteSmokeTests(unittest.TestCase):
    def test_pipe_only_case_converges(self) -> None:
        result = solve(
            build_steady_water_network_no_fittings_case(),
            SolverSettings(
                turbulent_iterations=60,
                pressure_relaxation=1.0,
                colebrook_residual_tolerance=1e-4,
                pressure_correction_abs_tolerance_pa=1e-3,
                pressure_correction_rel_tolerance=1e-8,
            ),
        )
        self.assertTrue(result.converged)

    def test_fittings_case_converges(self) -> None:
        result = solve(
            build_steady_water_network_case(),
            SolverSettings(
                turbulent_iterations=60,
                pressure_relaxation=1.0,
                colebrook_residual_tolerance=1e-4,
                pressure_correction_abs_tolerance_pa=1e-3,
                pressure_correction_rel_tolerance=1e-8,
            ),
        )
        self.assertTrue(result.converged)

    def test_elevation_case_converges(self) -> None:
        result = solve(
            build_steady_water_network_aggressive_elevation_case(),
            SolverSettings(
                turbulent_iterations=60,
                pressure_relaxation=1.0,
                colebrook_residual_tolerance=1e-4,
                pressure_correction_abs_tolerance_pa=1e-3,
                pressure_correction_rel_tolerance=1e-8,
            ),
        )
        self.assertTrue(result.converged)

    def test_inlet_flow_case_converges(self) -> None:
        result = solve(
            build_steady_water_network_inlet_flow_boundary_case(),
            SolverSettings(
                turbulent_iterations=60,
                pressure_relaxation=1.0,
                colebrook_residual_tolerance=1e-4,
                pressure_correction_abs_tolerance_pa=1e-3,
                pressure_correction_rel_tolerance=1e-6,
            ),
        )
        self.assertTrue(result.converged)

    def test_outlet_flow_case_converges(self) -> None:
        result = solve(
            build_steady_water_network_aggressive_elevation_outlet_flow_case(),
            SolverSettings(
                laminar_iterations=20,
                turbulent_iterations=120,
                pressure_relaxation=0.3,
                colebrook_residual_tolerance=1e-4,
                pressure_correction_abs_tolerance_pa=1e-3,
                pressure_correction_rel_tolerance=1e-6,
            ),
        )
        self.assertTrue(result.converged)

    def test_two_flow_boundaries_case_converges(self) -> None:
        result = solve(
            build_steady_water_network_two_flow_boundaries_case(),
            SolverSettings(
                turbulent_iterations=60,
                pressure_relaxation=1.0,
                colebrook_residual_tolerance=1e-4,
                pressure_correction_abs_tolerance_pa=1e-3,
                pressure_correction_rel_tolerance=1e-7,
            ),
        )
        self.assertTrue(result.converged)

    def test_three_reservoir_junction_converges(self) -> None:
        solver = SteadyIsothermalIncompressibleSolver(
            settings=SolverSettings(
                turbulent_iterations=300,
                pressure_relaxation=0.5,
                pressure_correction_abs_tolerance_pa=1e-3,
                pressure_correction_rel_tolerance=1e-8,
            ),
            turbulent_pipe_correlation=ColebrookPipeCorrelation(),
        )
        result = solver.solve(build_three_reservoir_junction_case())
        self.assertTrue(result.converged)
        junction_pa = result.node_pressures_pa[2]
        junction_head_m = junction_pa / (998.25 * 9.81)
        self.assertAlmostEqual(junction_head_m, 25.26, delta=0.01)
        flows = {cf.label: cf.volumetric_flow_m3_per_h for cf in result.component_flows}
        self.assertAlmostEqual(flows["Pipe 1->2"], 112.35, delta=0.1)
        self.assertAlmostEqual(flows["Pipe 2->3"], 62.36, delta=0.1)
        self.assertAlmostEqual(flows["Pipe 2->4"], 49.99, delta=0.1)
