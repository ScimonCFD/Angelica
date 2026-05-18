from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path
from types import SimpleNamespace

SRC_ROOT = Path(__file__).resolve().parents[1] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.closures.friction import ColebrookPipeCorrelation


class ColebrookFrictionTests(unittest.TestCase):
    def test_solve_colebrook_returns_positive_factor_with_transformed_strategy(self) -> None:
        pipe_state = SimpleNamespace(
            reynolds=2.5e5,
            component=SimpleNamespace(
                absolute_roughness_m=4.5e-5,
                diameter_m=0.3048,
            ),
        )
        correlation = ColebrookPipeCorrelation()

        friction_factor = correlation.solve_colebrook(
            pipe_state,
            initial_guess=1e-6,
            tolerance=1e-8,
            colebrook_friction_strategy="transformed",
            friction_factor_method="newton",
            friction_factor_max_iterations=80,
        )

        self.assertTrue(math.isfinite(friction_factor))
        self.assertGreater(friction_factor, 0.0)

    def test_solve_colebrook_returns_positive_factor_with_direct_strategy(self) -> None:
        pipe_state = SimpleNamespace(
            reynolds=2.5e5,
            component=SimpleNamespace(
                absolute_roughness_m=4.5e-5,
                diameter_m=0.3048,
            ),
        )
        correlation = ColebrookPipeCorrelation()

        friction_factor = correlation.solve_colebrook(
            pipe_state,
            initial_guess=1e-6,
            tolerance=1e-8,
            colebrook_friction_strategy="direct",
            friction_factor_method="newton",
            friction_factor_max_iterations=80,
        )

        self.assertTrue(math.isfinite(friction_factor))
        self.assertGreater(friction_factor, 0.0)

    def test_evaluate_colebrook_handles_non_positive_trial_value(self) -> None:
        pipe_state = SimpleNamespace(
            reynolds=5.0e4,
            component=SimpleNamespace(
                absolute_roughness_m=4.5e-5,
                diameter_m=0.1,
            ),
        )
        correlation = ColebrookPipeCorrelation()

        residual = correlation.evaluate_colebrook(pipe_state, -1.0e-3)

        self.assertTrue(math.isfinite(residual))


if __name__ == "__main__":
    unittest.main()
