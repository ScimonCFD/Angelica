from __future__ import annotations

import sys
import unittest
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[1] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.properties.dead_oil import dead_oil_density_kg_per_m3, dead_oil_viscosity_pa_s


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
