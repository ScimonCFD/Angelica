from __future__ import annotations

import sys
import unittest
from dataclasses import dataclass
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[1] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica import (
    CompressibleFluid,
    IdealGasEOS,
    NetworkCase,
    Pipe,
    PengRobinsonEOS,
    PressureBoundary,
    SteadyCompressibleSolver,
    lee_gonzalez_eakin_viscosity,
)
from angelica.solvers import CompressibleSolverSettings

# Methane critical properties (NIST)
_TC_CH4 = 190.564       # K
_PC_CH4 = 4_599_200.0   # Pa
_W_CH4  = 0.01141       # acentric factor
_M_CH4  = 0.016043      # kg/mol

_R = 8.314          # J/(mol·K)
_M_AIR = 0.028964   # kg/mol — dry air
_MU_AIR = 1.81e-5   # Pa·s  — air at 20 °C


# ── helpers ──────────────────────────────────────────────────────────────────

def _air_eos():
    return IdealGasEOS(molecular_weight_kg_per_mol=_M_AIR)


def _air_fluid(reference_pressure_pa=101_325.0):
    return CompressibleFluid.from_constants(
        eos=_air_eos(),
        viscosity_pa_s=_MU_AIR,
        specific_heat_j_per_kg_k=1005.0,
        thermal_conductivity_w_per_m_k=0.0257,
        reference_pressure_pa=reference_pressure_pa,
    )


def _single_pipe_case(p_in_pa, p_out_pa, diameter=0.05, length=100.0):
    """Source (node 1) → Pipe → Sink (node 2)."""
    return NetworkCase(
        name="compressible single pipe",
        fluid_model=_air_fluid(reference_pressure_pa=p_out_pa),
        pressure_inlets=(PressureBoundary(node_id=1, pressure_pa=p_in_pa),),
        pressure_outlets=(PressureBoundary(node_id=2, pressure_pa=p_out_pa),),
        components=(
            Pipe(
                start_node=1,
                end_node=2,
                diameter_m=diameter,
                length_m=length,
                absolute_roughness_m=0.000045,
            ),
        ),
        node_ids=(1, 2),
    )


def _two_pipe_case(p_in_pa, p_out_pa):
    """Source (1) → Pipe → Junction (2) → Pipe → Sink (3)."""
    return NetworkCase(
        name="compressible two pipes in series",
        fluid_model=_air_fluid(reference_pressure_pa=p_out_pa),
        pressure_inlets=(PressureBoundary(node_id=1, pressure_pa=p_in_pa),),
        pressure_outlets=(PressureBoundary(node_id=3, pressure_pa=p_out_pa),),
        components=(
            Pipe(start_node=1, end_node=2, diameter_m=0.05, length_m=50.0, absolute_roughness_m=0.000045),
            Pipe(start_node=2, end_node=3, diameter_m=0.05, length_m=50.0, absolute_roughness_m=0.000045),
        ),
        node_ids=(1, 2, 3),
    )


@dataclass
class _FakeNode:
    pressure_pa: float | None


@dataclass
class _FakeLinkState:
    start_node: _FakeNode
    end_node: _FakeNode
    temperature_c: float | None = None


def _solver():
    return SteadyCompressibleSolver(
        compressible_settings=CompressibleSolverSettings(
            max_density_iterations=30,
            density_rel_tolerance=1e-5,
        ),
    )


# ── IdealGasEOS ──────────────────────────────────────────────────────────────

class IdealGasEOSTests(unittest.TestCase):

    def test_density_at_standard_conditions(self):
        eos = _air_eos()
        p, t = 101_325.0, 20.0
        expected = p * _M_AIR / (_R * (t + 273.15))
        self.assertAlmostEqual(eos.density(p, t), expected, places=6)

    def test_density_proportional_to_pressure(self):
        eos = _air_eos()
        rho_1 = eos.density(100_000.0, 20.0)
        rho_2 = eos.density(200_000.0, 20.0)
        self.assertAlmostEqual(rho_2 / rho_1, 2.0, places=10)

    def test_density_inversely_proportional_to_absolute_temperature(self):
        eos = _air_eos()
        rho_cold = eos.density(101_325.0, 0.0)
        rho_hot  = eos.density(101_325.0, 100.0)
        t_cold_k = 273.15
        t_hot_k  = 373.15
        self.assertAlmostEqual(rho_cold / rho_hot, t_hot_k / t_cold_k, places=10)

    def test_invalid_molecular_weight_raises(self):
        with self.assertRaises(ValueError):
            IdealGasEOS(molecular_weight_kg_per_mol=0.0)
        with self.assertRaises(ValueError):
            IdealGasEOS(molecular_weight_kg_per_mol=-1.0)


# ── CompressibleFluid ─────────────────────────────────────────────────────────

class CompressibleFluidTests(unittest.TestCase):

    def test_density_uses_average_node_pressure(self):
        fluid = _air_fluid()
        p_start, p_end = 200_000.0, 100_000.0
        link = _FakeLinkState(
            start_node=_FakeNode(pressure_pa=p_start),
            end_node=_FakeNode(pressure_pa=p_end),
        )
        expected = _air_eos().density(0.5 * (p_start + p_end), 20.0)
        self.assertAlmostEqual(fluid.density_for_link(link), expected, places=10)

    def test_density_falls_back_to_reference_when_pressures_none(self):
        ref_p = 150_000.0
        fluid = _air_fluid(reference_pressure_pa=ref_p)
        link = _FakeLinkState(
            start_node=_FakeNode(pressure_pa=None),
            end_node=_FakeNode(pressure_pa=None),
        )
        expected = _air_eos().density(ref_p, 20.0)
        self.assertAlmostEqual(fluid.density_for_link(link), expected, places=10)

    def test_density_increases_with_node_pressure(self):
        fluid = _air_fluid()
        link_low  = _FakeLinkState(_FakeNode(100_000.0), _FakeNode(90_000.0))
        link_high = _FakeLinkState(_FakeNode(200_000.0), _FakeNode(190_000.0))
        self.assertGreater(
            fluid.density_for_link(link_high),
            fluid.density_for_link(link_low),
        )

    def test_viscosity_constant_when_using_from_constants(self):
        fluid = _air_fluid()
        link_low  = _FakeLinkState(_FakeNode(100_000.0), _FakeNode(90_000.0))
        link_high = _FakeLinkState(_FakeNode(500_000.0), _FakeNode(490_000.0))
        self.assertEqual(
            fluid.viscosity_for_link(link_low),
            fluid.viscosity_for_link(link_high),
        )


# ── PengRobinsonEOS ──────────────────────────────────────────────────────────

class PengRobinsonEOSTests(unittest.TestCase):

    def _methane_eos(self):
        return PengRobinsonEOS(
            molecular_weight_kg_per_mol=_M_CH4,
            critical_temperature_k=_TC_CH4,
            critical_pressure_pa=_PC_CH4,
            acentric_factor=_W_CH4,
        )

    def test_low_pressure_converges_to_ideal_gas(self):
        # At Pr ≈ 0.02 (near-atmospheric), PR and ideal gas agree to < 1 %
        eos_pr = self._methane_eos()
        eos_ig = IdealGasEOS(molecular_weight_kg_per_mol=_M_CH4)
        rho_pr = eos_pr.density(101_325.0, 20.0)
        rho_ig = eos_ig.density(101_325.0, 20.0)
        self.assertAlmostEqual(rho_pr / rho_ig, 1.0, delta=0.01)

    def test_density_increases_with_pressure(self):
        eos = self._methane_eos()
        rho_lo = eos.density(1_000_000.0, 20.0)
        rho_hi = eos.density(5_000_000.0, 20.0)
        self.assertGreater(rho_hi, rho_lo)

    def test_density_decreases_with_temperature(self):
        eos = self._methane_eos()
        rho_cold = eos.density(3_000_000.0, 10.0)
        rho_warm = eos.density(3_000_000.0, 60.0)
        self.assertGreater(rho_cold, rho_warm)

    def test_pr_denser_than_ideal_at_moderate_pressure(self):
        # For methane at Tr ≈ 1.5, Pr ≈ 0.65 (3 MPa, 20 °C), Z < 1
        # so PR density is higher than ideal gas density
        eos_pr = self._methane_eos()
        eos_ig = IdealGasEOS(molecular_weight_kg_per_mol=_M_CH4)
        rho_pr = eos_pr.density(3_000_000.0, 20.0)
        rho_ig = eos_ig.density(3_000_000.0, 20.0)
        self.assertGreater(rho_pr, rho_ig)

    def test_invalid_inputs_raise(self):
        with self.assertRaises(ValueError):
            PengRobinsonEOS(
                molecular_weight_kg_per_mol=-1.0,
                critical_temperature_k=_TC_CH4,
                critical_pressure_pa=_PC_CH4,
                acentric_factor=_W_CH4,
            )
        with self.assertRaises(ValueError):
            PengRobinsonEOS(
                molecular_weight_kg_per_mol=_M_CH4,
                critical_temperature_k=0.0,
                critical_pressure_pa=_PC_CH4,
                acentric_factor=_W_CH4,
            )


# ── lee_gonzalez_eakin_viscosity ──────────────────────────────────────────────

class LeeGonzalezEakinTests(unittest.TestCase):

    def _ch4_pr_eos(self):
        return PengRobinsonEOS(
            molecular_weight_kg_per_mol=_M_CH4,
            critical_temperature_k=_TC_CH4,
            critical_pressure_pa=_PC_CH4,
            acentric_factor=_W_CH4,
        )

    def test_atmospheric_value_close_to_known(self):
        # Methane at 1 atm, 20 °C: accepted value ≈ 1.1e-5 Pa·s
        eos = IdealGasEOS(molecular_weight_kg_per_mol=_M_CH4)
        mu_fn = lee_gonzalez_eakin_viscosity(eos, _M_CH4)
        mu = mu_fn(101_325.0, 20.0)
        self.assertAlmostEqual(mu, 1.1e-5, delta=0.1e-5)

    def test_viscosity_increases_with_temperature(self):
        # Gas viscosity increases with T (kinetic theory)
        eos = IdealGasEOS(molecular_weight_kg_per_mol=_M_CH4)
        mu_fn = lee_gonzalez_eakin_viscosity(eos, _M_CH4)
        self.assertGreater(mu_fn(101_325.0, 60.0), mu_fn(101_325.0, 20.0))

    def test_viscosity_increases_with_pressure(self):
        # Higher P → higher ρ → higher viscosity via LGE
        eos = self._ch4_pr_eos()
        mu_fn = lee_gonzalez_eakin_viscosity(eos, _M_CH4)
        self.assertGreater(mu_fn(5_000_000.0, 20.0), mu_fn(1_000_000.0, 20.0))

    def test_works_with_pr_eos(self):
        # PR EOS gives higher density than ideal → higher viscosity at 3 MPa
        eos_ig = IdealGasEOS(molecular_weight_kg_per_mol=_M_CH4)
        eos_pr = self._ch4_pr_eos()
        mu_ig = lee_gonzalez_eakin_viscosity(eos_ig, _M_CH4)(3_000_000.0, 20.0)
        mu_pr = lee_gonzalez_eakin_viscosity(eos_pr, _M_CH4)(3_000_000.0, 20.0)
        self.assertGreater(mu_pr, mu_ig)


# ── SteadyCompressibleSolver ──────────────────────────────────────────────────

class SteadyCompressibleSolverTests(unittest.TestCase):

    def test_single_pipe_converges(self):
        result = _solver().solve(_single_pipe_case(200_000.0, 100_000.0))
        self.assertTrue(result.converged)

    def test_zero_delta_p_gives_zero_flow(self):
        result = _solver().solve(_single_pipe_case(100_000.0, 100_000.0))
        mdot = result.component_flows[0].mass_flow_kg_per_s
        self.assertAlmostEqual(mdot, 0.0, delta=1e-6)

    def test_higher_delta_p_gives_higher_mass_flow(self):
        result_low  = _solver().solve(_single_pipe_case(110_000.0, 100_000.0))
        result_high = _solver().solve(_single_pipe_case(200_000.0, 100_000.0))
        mdot_low  = result_low.component_flows[0].mass_flow_kg_per_s
        mdot_high = result_high.component_flows[0].mass_flow_kg_per_s
        self.assertGreater(mdot_high, mdot_low)

    def test_mass_conservation_two_pipes_in_series(self):
        result = _solver().solve(_two_pipe_case(300_000.0, 100_000.0))
        self.assertTrue(result.converged)
        mdot_1 = result.component_flows[0].mass_flow_kg_per_s
        mdot_2 = result.component_flows[1].mass_flow_kg_per_s
        # Mass must be conserved at the junction node
        self.assertAlmostEqual(mdot_1, mdot_2, delta=abs(mdot_1) * 1e-4)

    def test_inlet_density_higher_than_outlet_density(self):
        eos = _air_eos()
        case = _single_pipe_case(200_000.0, 100_000.0)
        result = _solver().solve(case)
        self.assertTrue(result.converged)
        rho_in  = eos.density(result.node_pressures_pa[1], 20.0)
        rho_out = eos.density(result.node_pressures_pa[2], 20.0)
        self.assertGreater(rho_in, rho_out)

    def test_density_history_recorded(self):
        result = _solver().solve(_single_pipe_case(200_000.0, 100_000.0))
        self.assertGreater(len(result.density_history), 0)

    def test_density_history_converges(self):
        result = _solver().solve(_single_pipe_case(200_000.0, 100_000.0))
        self.assertLess(result.density_history[-1], 1e-4)


if __name__ == "__main__":
    unittest.main()
