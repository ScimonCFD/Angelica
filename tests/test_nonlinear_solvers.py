from __future__ import annotations

import sys
import unittest
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[1] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.numerics import (
    FixedPointSolver,
    NewtonSolver,
    NonlinearProblem,
    SecantSolver,
    build_nonlinear_solver,
)


class NonlinearSolverTests(unittest.TestCase):
    def test_fixed_point_solver_converges_on_simple_mapping(self) -> None:
        solver = FixedPointSolver()
        problem = NonlinearProblem(
            fixed_point_fn=lambda value: 0.5 * (value + 2.0 / value),
            residual_fn=lambda value: value * value - 2.0,
        )

        solution = solver.solve(problem, initial_value=1.0, tolerance=1e-8)

        self.assertAlmostEqual(solution, 2.0**0.5, places=6)

    def test_secant_solver_converges_on_simple_residual(self) -> None:
        solver = SecantSolver()
        problem = NonlinearProblem(
            residual_fn=lambda value: value * value - 2.0,
        )

        solution = solver.solve(problem, initial_value=1.0, tolerance=1e-8)

        self.assertAlmostEqual(solution, 2.0**0.5, places=6)

    def test_newton_solver_converges_on_simple_residual(self) -> None:
        solver = NewtonSolver()
        problem = NonlinearProblem(
            residual_fn=lambda value: value * value - 2.0,
            derivative_fn=lambda value: 2.0 * value,
        )

        solution = solver.solve(problem, initial_value=1.0, tolerance=1e-8)

        self.assertAlmostEqual(solution, 2.0**0.5, places=6)

    def test_factory_builds_known_nonlinear_solvers(self) -> None:
        self.assertIsInstance(build_nonlinear_solver("fixed_point"), FixedPointSolver)
        self.assertIsInstance(build_nonlinear_solver("secant"), SecantSolver)
        self.assertIsInstance(build_nonlinear_solver("newton"), NewtonSolver)

        with self.assertRaises(ValueError):
            build_nonlinear_solver("unknown")

    def test_factory_returns_same_instance_for_same_key(self) -> None:
        a = build_nonlinear_solver("newton", max_iterations=50)
        b = build_nonlinear_solver("newton", max_iterations=50)
        self.assertIs(a, b)

    def test_factory_returns_different_instances_for_different_max_iterations(self) -> None:
        a = build_nonlinear_solver("secant", max_iterations=30)
        b = build_nonlinear_solver("secant", max_iterations=80)
        self.assertIsNot(a, b)
        self.assertEqual(a.max_iterations, 30)
        self.assertEqual(b.max_iterations, 80)
