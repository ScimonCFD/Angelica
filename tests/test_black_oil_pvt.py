"""Unit tests for black-oil PVT correlations and BlackOilFluid model."""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[1] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.properties.black_oil import (
    BlackOilFluid,
    bubble_point_pa,
    gas_fvf,
    gas_viscosity_pa_s,
    live_oil_viscosity_pa_s,
    oil_fvf,
    solution_gor_m3_per_m3,
    water_fvf,
    water_viscosity_pa_s,
    z_factor_hall_yarborough,
    _M3M3_TO_SCFSTB,
    _PSIA_TO_PA,
)


class TestZFactorHallYarborough(unittest.TestCase):
    """Validate z-factor against published Standing-Katz chart values."""

    def test_z_near_unity_at_low_pressure(self):
        # At near-atmospheric conditions ideal gas behaviour: z ≈ 1
        z = z_factor_hall_yarborough(101_325.0, 20.0, gas_gravity=0.65)
        self.assertAlmostEqual(z, 1.0, delta=0.02)

    def test_z_decreases_then_increases_with_pressure(self):
        # Classic shape: z dips below 1 at moderate P, rises above at high P
        gas_gravity = 0.65
        T = 93.0  # °C ≈ 200 °F — typical reservoir temperature
        z_low  = z_factor_hall_yarborough(1e6,  T, gas_gravity)
        z_mid  = z_factor_hall_yarborough(7e6,  T, gas_gravity)
        z_high = z_factor_hall_yarborough(30e6, T, gas_gravity)
        self.assertGreater(z_low,  z_mid)   # z decreases initially
        self.assertGreater(z_high, z_mid)   # then increases

    def test_z_against_papay_reference(self):
        # Cross-validate vs Papay (1968) at T_pr=1.807, P_pr=2.984 (γ_g=0.65):
        #   z_Papay = 1 - 3.52*Ppr/10^(0.9813*Tpr) + 0.274*Ppr²/10^(0.8157*Tpr)
        #           ≈ 0.905
        # Hall-Yarborough is more accurate; agreement within 1–2 % is expected.
        z = z_factor_hall_yarborough(13.79e6, 93.3, gas_gravity=0.65)
        self.assertAlmostEqual(z, 0.90, delta=0.05)

    def test_z_positive(self):
        for P in [500e3, 2e6, 10e6, 30e6]:
            z = z_factor_hall_yarborough(P, 60.0, gas_gravity=0.70)
            self.assertGreater(z, 0.0)


class TestBubblePoint(unittest.TestCase):

    def test_zero_gor_gives_zero_pb(self):
        Pb = bubble_point_pa(0.0, 0.65, 35.0, 80.0)
        self.assertEqual(Pb, 0.0)

    def test_pb_increases_with_gor(self):
        # More gas → need higher pressure to keep it dissolved
        Pb_lo = bubble_point_pa(50.0,  0.65, 35.0, 80.0)
        Pb_hi = bubble_point_pa(200.0, 0.65, 35.0, 80.0)
        self.assertGreater(Pb_hi, Pb_lo)

    def test_pb_increases_with_temperature(self):
        # At higher T, gas is less soluble → Pb rises
        Pb_cold = bubble_point_pa(100.0, 0.65, 35.0, 40.0)
        Pb_hot  = bubble_point_pa(100.0, 0.65, 35.0, 100.0)
        self.assertGreater(Pb_hot, Pb_cold)

    def test_pb_standing_example(self):
        # Standing (1947): API=35, γ_g=0.65, T=150°F≈65.6°C, GOR=200 scf/STB
        gor_m3m3 = 200.0 / _M3M3_TO_SCFSTB
        Pb = bubble_point_pa(gor_m3m3, 0.65, 35.0, 65.6)
        Pb_psia = Pb / _PSIA_TO_PA
        # Expected ≈ 1030 psia; allow ±5 %
        self.assertAlmostEqual(Pb_psia, 1030.0, delta=60.0)


class TestSolutionGOR(unittest.TestCase):

    def test_rs_equals_gor_above_bubble_point(self):
        gor = 100.0  # m³/m³
        Pb  = bubble_point_pa(gor, 0.65, 35.0, 80.0)
        P_above = Pb * 1.5
        Rs = solution_gor_m3_per_m3(P_above, 80.0, 0.65, 35.0, gor)
        self.assertAlmostEqual(Rs, gor, delta=1e-6)

    def test_rs_less_than_gor_below_bubble_point(self):
        gor = 100.0
        Pb  = bubble_point_pa(gor, 0.65, 35.0, 80.0)
        P_below = Pb * 0.5
        Rs = solution_gor_m3_per_m3(P_below, 80.0, 0.65, 35.0, gor)
        self.assertLess(Rs, gor)

    def test_rs_increases_with_pressure(self):
        gor = 100.0
        Rs_lo = solution_gor_m3_per_m3(1e6, 80.0, 0.65, 35.0, gor)
        Rs_hi = solution_gor_m3_per_m3(5e6, 80.0, 0.65, 35.0, gor)
        self.assertGreater(Rs_hi, Rs_lo)

    def test_rs_zero_when_gor_zero(self):
        Rs = solution_gor_m3_per_m3(5e6, 80.0, 0.65, 35.0, 0.0)
        self.assertEqual(Rs, 0.0)


class TestOilFVF(unittest.TestCase):

    def test_bo_greater_than_one(self):
        # Bo always > 1: live oil is swollen relative to dead oil
        Rs = 50.0  # m³/m³
        Bo = oil_fvf(Rs, 80.0, 0.65, 35.0)
        self.assertGreater(Bo, 1.0)

    def test_bo_increases_with_rs(self):
        Bo_lo = oil_fvf(30.0, 80.0, 0.65, 35.0)
        Bo_hi = oil_fvf(150.0, 80.0, 0.65, 35.0)
        self.assertGreater(Bo_hi, Bo_lo)

    def test_bo_dead_oil_near_one(self):
        # Rs = 0 (dead oil): Bo should be close to 1
        Bo = oil_fvf(0.0, 20.0, 0.65, 35.0)
        self.assertAlmostEqual(Bo, 0.972, delta=0.05)


class TestGasFVF(unittest.TestCase):

    def test_bg_decreases_with_pressure(self):
        # At higher P, same mass of gas occupies less volume → Bg decreases
        Bg_lo = gas_fvf(1e6, 80.0, z=0.9)
        Bg_hi = gas_fvf(5e6, 80.0, z=0.85)
        self.assertGreater(Bg_lo, Bg_hi)

    def test_bg_ideal_gas_relation(self):
        # For z=1: Bg = P_sc * T / (T_sc * P)
        P  = 5e6
        T  = 80.0
        z  = 1.0
        Bg = gas_fvf(P, T, z)
        expected = 101_325.0 * (T + 273.15) / (288.706 * P)
        self.assertAlmostEqual(Bg, expected, delta=1e-10)


class TestWaterFVF(unittest.TestCase):

    def test_bw_near_one(self):
        Bw = water_fvf(5e6, 60.0)
        self.assertAlmostEqual(Bw, 1.0, delta=0.05)


class TestViscosities(unittest.TestCase):

    def test_live_oil_less_viscous_than_dead_oil(self):
        # Dissolved gas reduces oil viscosity
        from angelica.properties.dead_oil import dead_oil_viscosity_pa_s
        mu_dead = dead_oil_viscosity_pa_s(32.0, 60.0)
        mu_live = live_oil_viscosity_pa_s(mu_dead, rs_m3_per_m3=50.0)
        self.assertLess(mu_live, mu_dead)

    def test_gas_viscosity_increases_with_temperature_at_low_pressure(self):
        # Gas viscosity increases with T (unlike liquids)
        mu_cold = gas_viscosity_pa_s(101_325.0, 20.0, 0.65)
        mu_hot  = gas_viscosity_pa_s(101_325.0, 100.0, 0.65)
        self.assertGreater(mu_hot, mu_cold)

    def test_water_viscosity_decreases_with_temperature(self):
        mu_cold = water_viscosity_pa_s(20.0)
        mu_hot  = water_viscosity_pa_s(80.0)
        self.assertLess(mu_hot, mu_cold)

    def test_water_viscosity_at_20c(self):
        # Known value: ≈ 1 cP = 1e-3 Pa·s at 20 °C
        mu = water_viscosity_pa_s(20.0)
        self.assertAlmostEqual(mu * 1000.0, 1.0, delta=0.15)


class TestBlackOilFluid(unittest.TestCase):
    """Integration tests for the BlackOilFluid FluidModel."""

    def _make_fluid(self, gor=100.0, wor=0.3):
        return BlackOilFluid(
            api_gravity=35.0,
            gas_gravity=0.65,
            gor_sc_m3_per_m3=gor,
            wor_sc_m3_per_m3=wor,
            reference_pressure_pa=5e6,
            reference_temperature_c=60.0,
        )

    def test_density_positive(self):
        fluid = self._make_fluid()
        state = fluid.pvt(5e6, 60.0)
        self.assertGreater(state.mixture_density_kg_per_m3, 0.0)

    def test_viscosity_positive(self):
        fluid = self._make_fluid()
        state = fluid.pvt(5e6, 60.0)
        self.assertGreater(state.mixture_viscosity_pa_s, 0.0)

    def test_holdups_sum_to_one(self):
        fluid = self._make_fluid()
        # Two-phase case (P below Pb)
        Pb = bubble_point_pa(100.0, 0.65, 35.0, 60.0)
        state = fluid.pvt(Pb * 0.3, 60.0)
        total = state.holdup_oil + state.holdup_gas + state.holdup_water
        self.assertAlmostEqual(total, 1.0, places=10)

    def test_no_free_gas_above_bubble_point(self):
        fluid = self._make_fluid(gor=50.0, wor=0.0)
        Pb    = bubble_point_pa(50.0, 0.65, 35.0, 60.0)
        state = fluid.pvt(Pb * 2.0, 60.0)
        self.assertAlmostEqual(state.holdup_gas, 0.0, places=10)

    def test_free_gas_appears_below_bubble_point(self):
        fluid = self._make_fluid(gor=100.0, wor=0.0)
        Pb    = bubble_point_pa(100.0, 0.65, 35.0, 60.0)
        state = fluid.pvt(Pb * 0.5, 60.0)
        self.assertGreater(state.holdup_gas, 0.0)

    def test_dead_oil_has_no_gas(self):
        fluid = BlackOilFluid(
            api_gravity=35.0,
            gas_gravity=0.65,
            gor_sc_m3_per_m3=0.0,
            wor_sc_m3_per_m3=0.0,
        )
        state = fluid.pvt(5e6, 60.0)
        self.assertAlmostEqual(state.holdup_gas, 0.0, places=10)
        self.assertAlmostEqual(state.holdup_oil, 1.0, places=10)

    def test_density_for_link_uses_node_pressures(self):
        """density_for_link should read P from node pressures."""
        class FakeNode:
            pressure_pa = 5e6
        class FakeLink:
            start_node = FakeNode()
            end_node   = FakeNode()
            temperature_c = 60.0

        fluid = self._make_fluid()
        rho = fluid.density_for_link(FakeLink())
        self.assertGreater(rho, 0.0)

    def test_invalid_api_raises(self):
        with self.assertRaises(ValueError):
            BlackOilFluid(api_gravity=-200.0, gas_gravity=0.65,
                          gor_sc_m3_per_m3=100.0, wor_sc_m3_per_m3=0.0)

    def test_negative_gor_raises(self):
        with self.assertRaises(ValueError):
            BlackOilFluid(api_gravity=35.0, gas_gravity=0.65,
                          gor_sc_m3_per_m3=-1.0, wor_sc_m3_per_m3=0.0)


if __name__ == "__main__":
    unittest.main()
