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
)
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
