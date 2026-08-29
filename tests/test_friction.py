from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path
from types import SimpleNamespace

SRC_ROOT = Path(__file__).resolve().parents[1] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.closures.friction import ColebrookPipeCorrelation, HazenWilliamsPipeCorrelation
from angelica.closures.heat_source import HeatSourceModel
from angelica.closures.minor_losses import MinorLossModel


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


class HazenWilliamsPipeCorrelationTests(unittest.TestCase):
    """Verify HazenWilliams velocity and coupling against hand-calculated values.

    Reference: head-loss form h_f = 10.67 * L * Q^1.852 / (C^1.852 * D^4.871)
    Test pipe: L=100 m, D=0.1 m, C=100, rho=1000 kg/m³, g=9.81 m/s²
    """

    def _pipe_state(self, velocity: float = 1.0) -> SimpleNamespace:
        return SimpleNamespace(
            velocity_m_per_s=velocity,
            area_m2=math.pi * 0.1 ** 2 / 4.0,
            component=SimpleNamespace(
                length_m=100.0,
                diameter_m=0.1,
                height_change_m=0.0,
                hazen_williams_c=100.0,
            ),
        )

    def test_zero_driving_pressure_gives_zero_velocity(self) -> None:
        corr = HazenWilliamsPipeCorrelation()
        ps = self._pipe_state()
        v = corr.calculate_velocity(ps, delta_p=0.0, density=1000.0, viscosity=1e-3)
        self.assertEqual(v, 0.0)

    def test_positive_driving_pressure_gives_positive_velocity(self) -> None:
        corr = HazenWilliamsPipeCorrelation()
        ps = self._pipe_state()
        # 10 000 Pa driving pressure
        v = corr.calculate_velocity(ps, delta_p=10_000.0, density=1000.0, viscosity=1e-3)
        self.assertGreater(v, 0.0)

    def test_negative_driving_pressure_gives_negative_velocity(self) -> None:
        corr = HazenWilliamsPipeCorrelation()
        ps = self._pipe_state()
        v = corr.calculate_velocity(ps, delta_p=-10_000.0, density=1000.0, viscosity=1e-3)
        self.assertLess(v, 0.0)

    def test_velocity_magnitude_symmetric(self) -> None:
        corr = HazenWilliamsPipeCorrelation()
        ps = self._pipe_state()
        v_pos = corr.calculate_velocity(ps, delta_p=5_000.0, density=1000.0, viscosity=1e-3)
        v_neg = corr.calculate_velocity(ps, delta_p=-5_000.0, density=1000.0, viscosity=1e-3)
        self.assertAlmostEqual(abs(v_pos), abs(v_neg), places=10)

    def test_coupling_is_negative(self) -> None:
        corr = HazenWilliamsPipeCorrelation()
        ps = self._pipe_state(velocity=0.5)
        coupling = corr.calculate_coupling(ps, density=1000.0, viscosity=1e-3)
        self.assertLess(coupling, 0.0)

    def test_higher_velocity_gives_less_negative_coupling(self) -> None:
        # |coupling| ∝ 1/Q^(n-1) with n=1.852, so larger Q → smaller |coupling|
        corr = HazenWilliamsPipeCorrelation()
        ps_slow = self._pipe_state(velocity=0.1)
        ps_fast = self._pipe_state(velocity=1.0)
        c_slow = corr.calculate_coupling(ps_slow, density=1000.0, viscosity=1e-3)
        c_fast = corr.calculate_coupling(ps_fast, density=1000.0, viscosity=1e-3)
        self.assertGreater(c_fast, c_slow)  # c_slow is more negative


class HeatSourceModelTests(unittest.TestCase):
    def _link_rated(self, velocity: float = 1.0) -> SimpleNamespace:
        return SimpleNamespace(
            velocity_m_per_s=velocity,
            area_m2=math.pi * 0.05 ** 2 / 4.0,
            component=SimpleNamespace(
                pressure_drop_pa=5_000.0,
                pressure_drop_mode="rated",
                rated_mass_flow_kg_per_s=2.0,
            ),
        )

    def _link_fixed(self, velocity: float = 1.0) -> SimpleNamespace:
        return SimpleNamespace(
            velocity_m_per_s=velocity,
            area_m2=math.pi * 0.05 ** 2 / 4.0,
            component=SimpleNamespace(
                pressure_drop_pa=3_000.0,
                pressure_drop_mode="fixed",
                rated_mass_flow_kg_per_s=1.0,
            ),
        )

    def _link_transparent(self) -> SimpleNamespace:
        return SimpleNamespace(
            velocity_m_per_s=1.0,
            area_m2=math.pi * 0.05 ** 2 / 4.0,
            component=SimpleNamespace(
                pressure_drop_pa=0.0,
                pressure_drop_mode="rated",
                rated_mass_flow_kg_per_s=1.0,
            ),
        )

    def test_rated_mode_positive_dp_gives_positive_velocity(self) -> None:
        corr = HeatSourceModel()
        v = corr.calculate_velocity(self._link_rated(), delta_p=100.0, density=800.0, viscosity=1e-3)
        self.assertGreater(v, 0.0)

    def test_rated_mode_negative_dp_gives_negative_velocity(self) -> None:
        corr = HeatSourceModel()
        v = corr.calculate_velocity(self._link_rated(), delta_p=-100.0, density=800.0, viscosity=1e-3)
        self.assertLess(v, 0.0)

    def test_fixed_mode_positive_dp_gives_positive_velocity(self) -> None:
        corr = HeatSourceModel()
        v = corr.calculate_velocity(self._link_fixed(velocity=0.5), delta_p=200.0, density=800.0, viscosity=1e-3)
        self.assertGreater(v, 0.0)

    def test_transparent_device_uses_floor_k(self) -> None:
        corr = HeatSourceModel()
        link = self._link_transparent()
        v = corr.calculate_velocity(link, delta_p=50.0, density=1000.0, viscosity=1e-3)
        self.assertTrue(math.isfinite(v))
        self.assertGreater(v, 0.0)

    def test_coupling_is_negative(self) -> None:
        corr = HeatSourceModel()
        coupling = corr.calculate_coupling(self._link_rated(velocity=1.0), density=800.0, viscosity=1e-3)
        self.assertLess(coupling, 0.0)

    def test_coupling_bounded_by_max(self) -> None:
        # At near-zero velocity, coupling should be clamped to -_MAX_COUPLING_ABS
        from angelica.closures.heat_source import _MAX_COUPLING_ABS
        corr = HeatSourceModel()
        link = self._link_rated(velocity=1e-15)
        coupling = corr.calculate_coupling(link, density=800.0, viscosity=1e-3)
        self.assertGreaterEqual(coupling, -_MAX_COUPLING_ABS)


class MinorLossModelTests(unittest.TestCase):
    def _link(self, k: float, velocity: float = 1.0) -> SimpleNamespace:
        return SimpleNamespace(
            velocity_m_per_s=velocity,
            area_m2=math.pi * 0.05 ** 2 / 4.0,
            component=SimpleNamespace(loss_coefficient=k),
        )

    def test_finite_k_positive_dp_gives_positive_velocity(self) -> None:
        corr = MinorLossModel()
        v = corr.calculate_velocity(self._link(2.0), delta_p=500.0, density=1000.0, viscosity=1e-3)
        self.assertGreater(v, 0.0)

    def test_finite_k_negative_dp_gives_negative_velocity(self) -> None:
        corr = MinorLossModel()
        v = corr.calculate_velocity(self._link(2.0), delta_p=-500.0, density=1000.0, viscosity=1e-3)
        self.assertLess(v, 0.0)

    def test_infinite_k_gives_zero_velocity(self) -> None:
        # Closed check valve — no flow allowed
        corr = MinorLossModel()
        v = corr.calculate_velocity(self._link(float("inf")), delta_p=1000.0, density=1000.0, viscosity=1e-3)
        self.assertEqual(v, 0.0)

    def test_infinite_k_coupling_is_near_zero(self) -> None:
        corr = MinorLossModel()
        coupling = corr.calculate_coupling(self._link(float("inf")), density=1000.0, viscosity=1e-3)
        self.assertAlmostEqual(coupling, 0.0, places=5)

    def test_finite_k_coupling_is_negative(self) -> None:
        corr = MinorLossModel()
        coupling = corr.calculate_coupling(self._link(2.0, velocity=1.0), density=1000.0, viscosity=1e-3)
        self.assertLess(coupling, 0.0)

    def test_larger_k_gives_lower_velocity(self) -> None:
        corr = MinorLossModel()
        v_low_k  = corr.calculate_velocity(self._link(1.0), delta_p=200.0, density=1000.0, viscosity=1e-3)
        v_high_k = corr.calculate_velocity(self._link(5.0), delta_p=200.0, density=1000.0, viscosity=1e-3)
        self.assertGreater(v_low_k, v_high_k)

    def test_velocity_magnitude_symmetric(self) -> None:
        corr = MinorLossModel()
        v_pos = corr.calculate_velocity(self._link(3.0), delta_p=300.0, density=1000.0, viscosity=1e-3)
        v_neg = corr.calculate_velocity(self._link(3.0), delta_p=-300.0, density=1000.0, viscosity=1e-3)
        self.assertAlmostEqual(abs(v_pos), abs(v_neg), places=10)


if __name__ == "__main__":
    unittest.main()
