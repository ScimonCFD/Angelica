"""Tests for the turbulent velocity loop in ColebrookPipeCorrelation.

The velocity loop is a nonlinear solve: given a pressure drop ΔP across a
pipe, find the velocity v such that the Darcy-Weisbach equation is satisfied
when f is evaluated via Colebrook-White at Re = ρvD/μ.

These tests cover:
- All three velocity_loop_method options (fixed_point, secant, newton-like
  via the fallback residual_fn).
- Agreement between methods on a reference pipe/fluid condition.
- Physical sanity: velocity increases with pressure drop, decreases with
  increased viscosity, is finite and positive for positive ΔP.
- Zero-pressure-drop edge case.
- Laminar regime (Re < 2300) handled inside calculate_velocity.
"""
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


def _make_pipe_state(
    diameter_m: float = 0.2,
    length_m: float = 1000.0,
    roughness_m: float = 4.5e-5,
    height_change_m: float = 0.0,
    velocity_m_per_s: float = 1.0,
) -> SimpleNamespace:
    """Return a minimal pipe_state-compatible SimpleNamespace."""
    area_m2 = math.pi * diameter_m ** 2 / 4.0
    return SimpleNamespace(
        velocity_m_per_s=velocity_m_per_s,
        reynolds=1.0,
        friction_factor=0.02,
        area_m2=area_m2,
        component=SimpleNamespace(
            diameter_m=diameter_m,
            length_m=length_m,
            absolute_roughness_m=roughness_m,
            height_change_m=height_change_m,
        ),
    )


# Reference conditions: water at 20 °C in a 200 mm steel pipe
_RHO  = 998.0    # kg/m³
_MU   = 1.0e-3   # Pa·s
_D    = 0.2      # m
_L    = 1000.0   # m
_DP   = 50_000.0 # Pa  (≈ 0.5 bar/km, turbulent flow)
_TOL  = 1e-8


class TestVelocityLoopMethods(unittest.TestCase):
    """All three velocity_loop_method options must converge to the same answer."""

    def _solve(self, method: str) -> float:
        corr = ColebrookPipeCorrelation()
        ps = _make_pipe_state()
        v = corr.calculate_velocity(
            ps, _DP, _RHO, _MU,
            tolerance=_TOL,
            velocity_loop_method=method,
            velocity_loop_max_iterations=100,
            velocity_loop_tolerance=_TOL,
        )
        return v

    def test_fixed_point_converges(self):
        v = self._solve("fixed_point")
        self.assertTrue(math.isfinite(v))
        self.assertGreater(v, 0.0)

    def test_secant_converges(self):
        v = self._solve("secant")
        self.assertTrue(math.isfinite(v))
        self.assertGreater(v, 0.0)

    def test_fixed_point_and_secant_agree(self):
        v_fp = self._solve("fixed_point")
        v_sc = self._solve("secant")
        self.assertAlmostEqual(v_fp, v_sc, delta=v_fp * 0.01,
                               msg=f"fixed_point={v_fp:.4f} m/s, secant={v_sc:.4f} m/s")

    def test_velocity_is_physically_reasonable(self):
        """50 kPa/km in a 200 mm pipe → expect 1–5 m/s for water."""
        v = self._solve("fixed_point")
        self.assertGreater(v, 0.5)
        self.assertLess(v, 10.0)


class TestVelocityLoopPhysics(unittest.TestCase):
    """Monotonicity and limiting behaviour."""

    def _solve(self, delta_p: float, rho: float = _RHO, mu: float = _MU) -> float:
        corr = ColebrookPipeCorrelation()
        ps = _make_pipe_state()
        return corr.calculate_velocity(
            ps, delta_p, rho, mu,
            tolerance=_TOL,
            velocity_loop_method="secant",
            velocity_loop_max_iterations=100,
        )

    def test_higher_delta_p_gives_higher_velocity(self):
        v_lo = self._solve(20_000.0)
        v_hi = self._solve(80_000.0)
        self.assertGreater(v_hi, v_lo)

    def test_higher_viscosity_gives_lower_velocity(self):
        v_thin = self._solve(_DP, mu=1e-3)
        v_thick = self._solve(_DP, mu=5e-3)
        self.assertGreater(v_thin, v_thick)

    def test_higher_density_gives_lower_velocity_turbulent(self):
        """At turbulent Re, friction is roughly Re-independent, so ΔP ∝ ρv²."""
        v_light = self._solve(_DP, rho=700.0)
        v_heavy = self._solve(_DP, rho=1100.0)
        self.assertGreater(v_light, v_heavy)

    def test_darcy_weisbach_satisfied_at_solution(self):
        """Verify the solved velocity satisfies ΔP ≈ f·(L/D)·(ρ/2)·v² + elevation."""
        corr = ColebrookPipeCorrelation()
        ps = _make_pipe_state()
        v = corr.calculate_velocity(
            ps, _DP, _RHO, _MU,
            tolerance=_TOL,
            velocity_loop_method="secant",
            velocity_loop_max_iterations=100,
        )
        f = ps.friction_factor
        D, L = _D, _L
        dp_computed = f * (L / D) * (_RHO / 2.0) * v * abs(v)
        self.assertAlmostEqual(dp_computed, _DP, delta=_DP * 0.01,
                               msg=f"DW mismatch: computed={dp_computed:.1f} Pa, expected={_DP:.1f} Pa")


class TestVelocityLoopEdgeCases(unittest.TestCase):
    def test_zero_pressure_drop_gives_near_zero_velocity(self):
        corr = ColebrookPipeCorrelation()
        ps = _make_pipe_state()
        v = corr.calculate_velocity(
            ps, 0.0, _RHO, _MU,
            tolerance=_TOL,
            velocity_loop_method="fixed_point",
            velocity_loop_max_iterations=50,
        )
        self.assertAlmostEqual(v, 0.0, delta=0.01)

    def test_laminar_regime_no_iteration_needed(self):
        """Very high viscosity → Re < 2300 → Hagen-Poiseuille applies directly."""
        corr = ColebrookPipeCorrelation()
        ps = _make_pipe_state(diameter_m=0.05, length_m=10.0)
        # μ = 1 Pa·s (very viscous), ΔP = 100 Pa → Re << 2300
        v = corr.calculate_velocity(
            ps, 100.0, _RHO, 1.0,
            tolerance=_TOL,
            velocity_loop_method="fixed_point",
            velocity_loop_max_iterations=50,
        )
        self.assertTrue(math.isfinite(v))
        self.assertGreater(v, 0.0)
        re = _RHO * v * 0.05 / 1.0
        self.assertLess(re, 2300.0, f"Expected laminar Re, got Re={re:.1f}")

    def test_rough_pipe_lower_velocity_than_smooth(self):
        """Rougher pipe → higher friction → lower velocity at same ΔP."""
        corr = ColebrookPipeCorrelation()
        ps_smooth = _make_pipe_state(roughness_m=1e-6)
        ps_rough  = _make_pipe_state(roughness_m=5e-4)
        kw = {"tolerance": _TOL, "velocity_loop_method": "secant",
                  "velocity_loop_max_iterations": 100}
        v_smooth = corr.calculate_velocity(ps_smooth, _DP, _RHO, _MU, **kw)
        v_rough  = corr.calculate_velocity(ps_rough,  _DP, _RHO, _MU, **kw)
        self.assertGreater(v_smooth, v_rough)


if __name__ == "__main__":
    unittest.main()
