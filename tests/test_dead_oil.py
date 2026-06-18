from __future__ import annotations

import sys
import unittest
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[1] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.properties.dead_oil import (
    build_thermal_dead_oil,
    dead_oil_density_kg_per_m3,
    dead_oil_specific_heat_j_per_kg_k,
    dead_oil_thermal_conductivity_w_per_m_k,
    dead_oil_viscosity_pa_s,
)


class DeadOilDensityTests(unittest.TestCase):
    def test_30api_density(self) -> None:
        # SG = 141.5 / (30 + 131.5) = 0.8762 → ρ = 0.8762 × 999.064
        rho = dead_oil_density_kg_per_m3(30.0)
        self.assertAlmostEqual(rho, 875.5, delta=0.5)

    def test_40api_density(self) -> None:
        # SG = 141.5 / 171.5 = 0.8250 → ρ = 0.8250 × 999.064
        rho = dead_oil_density_kg_per_m3(40.0)
        self.assertAlmostEqual(rho, 824.5, delta=0.5)

    def test_density_decreases_with_api(self) -> None:
        self.assertGreater(
            dead_oil_density_kg_per_m3(20.0),
            dead_oil_density_kg_per_m3(45.0),
        )

    def test_10api_heavy_crude(self) -> None:
        rho = dead_oil_density_kg_per_m3(10.0)
        self.assertGreater(rho, 999.0)

    def test_invalid_api_raises(self) -> None:
        with self.assertRaises(ValueError):
            dead_oil_density_kg_per_m3(-200.0)


class DeadOilViscosityTests(unittest.TestCase):
    def test_viscosity_is_positive(self) -> None:
        self.assertGreater(dead_oil_viscosity_pa_s(35.0, 60.0), 0.0)

    def test_viscosity_decreases_with_temperature(self) -> None:
        mu_cold = dead_oil_viscosity_pa_s(30.0, 20.0)
        mu_hot = dead_oil_viscosity_pa_s(30.0, 80.0)
        self.assertGreater(mu_cold, mu_hot)

    def test_viscosity_decreases_with_api(self) -> None:
        mu_heavy = dead_oil_viscosity_pa_s(20.0, 60.0)
        mu_light = dead_oil_viscosity_pa_s(45.0, 60.0)
        self.assertGreater(mu_heavy, mu_light)

    def test_beggs_robinson_30api_60f(self) -> None:
        # At 60°F (15.56°C), 30°API: z = 3.0324 - 0.02023*30 - 1.163*log10(60)
        # Reference: ~180 cP
        mu_pa_s = dead_oil_viscosity_pa_s(30.0, 15.56)
        mu_cp = mu_pa_s * 1000.0
        self.assertGreater(mu_cp, 50.0)
        self.assertLess(mu_cp, 500.0)

    def test_light_crude_low_viscosity(self) -> None:
        # 45°API crude at 60°C — should be well under 10 cP
        mu_pa_s = dead_oil_viscosity_pa_s(45.0, 60.0)
        self.assertLess(mu_pa_s * 1000.0, 10.0)

    def test_invalid_temperature_raises(self) -> None:
        with self.assertRaises(ValueError):
            dead_oil_viscosity_pa_s(30.0, -300.0)


class DeadOilSpecificHeatTests(unittest.TestCase):
    def test_specific_heat_is_positive(self) -> None:
        self.assertGreater(dead_oil_specific_heat_j_per_kg_k(30.0, 60.0), 0.0)

    def test_specific_heat_increases_with_temperature(self) -> None:
        cp_cold = dead_oil_specific_heat_j_per_kg_k(30.0, 20.0)
        cp_hot = dead_oil_specific_heat_j_per_kg_k(30.0, 100.0)
        self.assertGreater(cp_hot, cp_cold)

    def test_specific_heat_increases_with_api(self) -> None:
        # Lighter crude (higher API) has higher cp
        cp_heavy = dead_oil_specific_heat_j_per_kg_k(20.0, 60.0)
        cp_light = dead_oil_specific_heat_j_per_kg_k(45.0, 60.0)
        self.assertGreater(cp_light, cp_heavy)

    def test_30api_at_60c_in_expected_range(self) -> None:
        # Watson-Nelson: typical crude oil cp is 1600–2500 J/(kg·K)
        cp = dead_oil_specific_heat_j_per_kg_k(30.0, 60.0)
        self.assertGreater(cp, 1600.0)
        self.assertLess(cp, 2500.0)


class DeadOilThermalConductivityTests(unittest.TestCase):
    def test_thermal_conductivity_is_positive(self) -> None:
        self.assertGreater(dead_oil_thermal_conductivity_w_per_m_k(30.0, 60.0), 0.0)

    def test_thermal_conductivity_decreases_with_temperature(self) -> None:
        k_cold = dead_oil_thermal_conductivity_w_per_m_k(30.0, 20.0)
        k_hot = dead_oil_thermal_conductivity_w_per_m_k(30.0, 100.0)
        self.assertGreater(k_cold, k_hot)

    def test_30api_at_20c_in_expected_range(self) -> None:
        # Cragoe (1929): typical crude oil k is 0.10–0.16 W/(m·K)
        k = dead_oil_thermal_conductivity_w_per_m_k(30.0, 20.0)
        self.assertGreater(k, 0.08)
        self.assertLess(k, 0.20)


class BuildThermalDeadOilTests(unittest.TestCase):
    def test_returns_thermal_fluid(self) -> None:
        from angelica.properties.thermal_fluid import ThermalFluid
        fluid = build_thermal_dead_oil(32.0)
        self.assertIsInstance(fluid, ThermalFluid)

    def test_density_matches_standard_correlation(self) -> None:
        fluid = build_thermal_dead_oil(32.0)
        # density is constant w.r.t. T (incompressible)
        rho_ref = dead_oil_density_kg_per_m3(32.0)

        class _FakeLink:
            temperature_c = 40.0

        self.assertAlmostEqual(fluid.density_for_link(_FakeLink()), rho_ref, places=3)

    def test_viscosity_decreases_with_temperature(self) -> None:
        fluid = build_thermal_dead_oil(32.0)

        class _Link:
            def __init__(self, t): self.temperature_c = t

        self.assertGreater(
            fluid.viscosity_for_link(_Link(20.0)),
            fluid.viscosity_for_link(_Link(80.0)),
        )

    def test_specific_heat_and_conductivity_accessible(self) -> None:
        fluid = build_thermal_dead_oil(32.0)

        class _Link:
            temperature_c = 60.0

        self.assertGreater(fluid.specific_heat_for_link(_Link()), 0.0)
        self.assertGreater(fluid.thermal_conductivity_for_link(_Link()), 0.0)
