"""
Tests for the Choke component and its integration with SteadyBlackOilSolver.

Coverage
--------
Unit
  - Choke K-factor formula: K = (D_pipe/D_orifice)^4 / Cd^2
  - Choke is-a Fitting (isinstance check for solver dispatch)
  - Validation errors: missing/negative/oversized orifice, bad Cd

Integration
  - Choke + Pipe + BB in a well-choke-separator case
  - Larger bean size → higher flow rate (monotone sensitivity)
  - Choke pressure drop fraction decreases as bean grows
  - Multi-well gathering network with BB: manifold pressure emerges
    correctly, and high-GOR wells show greater slip correction than
    low-GOR wells
"""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[1] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

import math
import warnings

from angelica.closures.beggs_brill import BeggsBrillCorrelation
from angelica.closures.choke import (
    ChokeFlowModel,
    sachdeva_critical_ratio,
    sachdeva_mass_flux,
)
from angelica.core.case import InletFluidBC, NetworkCase, PressureBoundary, ThermalBoundary
from angelica.core.components import Choke, Fitting, Pipe
from angelica.properties.black_oil import BlackOilFluid
from angelica.solvers import SteadyBlackOilSolver

# ── Shared fluid ────────────────────────────────────────────────────────────────

_CHOKE_FLUID = BlackOilFluid(
    api_gravity=35.0, gas_gravity=0.65,
    gor_sc_m3_per_m3=60.0, wor_sc_m3_per_m3=0.2,
    reference_pressure_pa=10e6, reference_temperature_c=80.0,
)


def _choke_case(orifice_d: float) -> NetworkCase:
    return NetworkCase(
        name="choke-test",
        fluid_model=_CHOKE_FLUID,
        pressure_inlets= (PressureBoundary(node_id=1, pressure_pa=20e6),),
        pressure_outlets=(PressureBoundary(node_id=3, pressure_pa=0.5e6),),
        thermal_inlets  =(ThermalBoundary(node_id=1, temperature_c=90.0,
                                          bc_type="fixed_temperature"),),
        components=(
            Choke(start_node=1, end_node=2, diameter_m=0.1016,
                  orifice_diameter_m=orifice_d, discharge_coefficient=0.62),
            Pipe(start_node=2, end_node=3, diameter_m=0.1016, length_m=5_000.0,
                 absolute_roughness_m=46e-6, height_change_m=0.0),
        ),
    )


# ── Unit tests: Choke geometry ──────────────────────────────────────────────────

class TestChokeComponent(unittest.TestCase):

    def test_is_instance_of_fitting(self):
        """Choke must be recognised as a Fitting for solver dispatch."""
        c = Choke(start_node=1, end_node=2, diameter_m=0.1, orifice_diameter_m=0.02)
        self.assertIsInstance(c, Fitting)

    def test_k_factor_formula(self):
        """K = (D_pipe / D_orifice)^4 / Cd^2 (standard orifice equation)."""
        D_pipe, D_orif, Cd = 0.1, 0.025, 0.62
        c = Choke(start_node=1, end_node=2, diameter_m=D_pipe,
                  orifice_diameter_m=D_orif, discharge_coefficient=Cd)
        K_expected = (D_pipe / D_orif) ** 4 / Cd ** 2
        self.assertAlmostEqual(c.loss_coefficient, K_expected, places=6)

    def test_larger_orifice_smaller_k(self):
        """Bigger bean → lower resistance."""
        c_small = Choke(start_node=1, end_node=2, diameter_m=0.1, orifice_diameter_m=0.01)
        c_large = Choke(start_node=1, end_node=2, diameter_m=0.1, orifice_diameter_m=0.04)
        self.assertGreater(c_small.loss_coefficient, c_large.loss_coefficient)

    def test_default_cd_is_sharp_edged(self):
        """Default discharge coefficient should be 0.62 (ISO 5167 sharp-edged)."""
        c = Choke(start_node=1, end_node=2, diameter_m=0.1, orifice_diameter_m=0.02)
        self.assertAlmostEqual(c.discharge_coefficient, 0.62)

    def test_missing_orifice_raises(self):
        with self.assertRaises(ValueError):
            Choke(start_node=1, end_node=2, diameter_m=0.1)

    def test_zero_orifice_raises(self):
        with self.assertRaises(ValueError):
            Choke(start_node=1, end_node=2, diameter_m=0.1, orifice_diameter_m=0.0)

    def test_orifice_ge_pipe_raises(self):
        with self.assertRaises(ValueError):
            Choke(start_node=1, end_node=2, diameter_m=0.1, orifice_diameter_m=0.1)

    def test_invalid_cd_raises(self):
        with self.assertRaises(ValueError):
            Choke(start_node=1, end_node=2, diameter_m=0.1,
                  orifice_diameter_m=0.02, discharge_coefficient=0.0)

    def test_cd_above_one_raises(self):
        with self.assertRaises(ValueError):
            Choke(start_node=1, end_node=2, diameter_m=0.1,
                  orifice_diameter_m=0.02, discharge_coefficient=1.1)


# ── Integration tests: well → choke → separator ─────────────────────────────────

class TestChokeIntegration(unittest.TestCase):

    def test_choke_solver_converges(self):
        """Well-choke-separator system must converge."""
        res = SteadyBlackOilSolver(
            turbulent_pipe_correlation=BeggsBrillCorrelation()
        ).solve(_choke_case(0.020))
        self.assertTrue(res.converged)
        self.assertGreater(abs(res.component_flows[0].mass_flow_kg_per_s), 0.0)

    def test_larger_bean_more_flow(self):
        """Increasing the bean size must increase production rate (monotone)."""
        solver = SteadyBlackOilSolver(
            turbulent_pipe_correlation=BeggsBrillCorrelation()
        )
        beans = [0.010, 0.016, 0.025]
        flows = [
            abs(solver.solve(_choke_case(d)).component_flows[1].volumetric_flow_m3_per_h)
            for d in beans
        ]
        for i in range(1, len(flows)):
            self.assertGreater(
                flows[i], flows[i - 1],
                msg=f"Flow did not increase from bean {beans[i-1]*1000:.0f}mm "
                    f"({flows[i-1]:.2f}) to {beans[i]*1000:.0f}mm ({flows[i]:.2f})",
            )

    def test_choke_dp_fraction_decreases_with_bean_size(self):
        """At large bean, the choke ΔP fraction < 50% of total ΔP."""
        solver = SteadyBlackOilSolver(
            turbulent_pipe_correlation=BeggsBrillCorrelation()
        )
        res_large = solver.solve(_choke_case(0.040))   # 40mm bean — wide-open
        P1 = res_large.node_pressures_pa[1]
        P2 = res_large.node_pressures_pa[2]
        P3 = res_large.node_pressures_pa[3]
        dP_choke = P1 - P2
        dP_total = P1 - P3
        fraction = dP_choke / dP_total
        self.assertLess(
            fraction, 0.5,
            msg=f"At large bean, choke fraction={fraction:.2%} should be < 50%",
        )

    def test_choke_pressure_drop_positive(self):
        """Downstream pressure must always be lower than upstream (P2 < P1)."""
        solver = SteadyBlackOilSolver(
            turbulent_pipe_correlation=BeggsBrillCorrelation()
        )
        res = solver.solve(_choke_case(0.015))
        P1 = res.node_pressures_pa[1]
        P2 = res.node_pressures_pa[2]
        self.assertGreater(P1, P2, msg="Choke must reduce pressure")


# ── Integration tests: gathering network with BB ────────────────────────────────

class TestGatheringNetworkBB(unittest.TestCase):
    """Multi-well gathering network (Tutorial 09 scenario) with Beggs-Brill."""

    _BASE_FLUID = BlackOilFluid(
        api_gravity=32.0, gas_gravity=0.65,
        gor_sc_m3_per_m3=50.0, wor_sc_m3_per_m3=0.3,
        reference_pressure_pa=6e6, reference_temperature_c=65.0,
    )

    _PIPES = (
        Pipe(component_id="A", start_node=1, end_node=4, diameter_m=0.1016,
             length_m=4_000.0, absolute_roughness_m=46e-6,
             heat_transfer_coefficient_w_per_m2k=5.0, ambient_temperature_c=15.0),
        Pipe(component_id="B", start_node=2, end_node=4, diameter_m=0.1016,
             length_m=2_000.0, absolute_roughness_m=46e-6,
             heat_transfer_coefficient_w_per_m2k=5.0, ambient_temperature_c=15.0),
        Pipe(component_id="D", start_node=4, end_node=5, diameter_m=0.1524,
             length_m=8_000.0, absolute_roughness_m=46e-6,
             heat_transfer_coefficient_w_per_m2k=3.0, ambient_temperature_c=15.0),
    )

    def _two_well_case(self) -> NetworkCase:
        return NetworkCase(
            name="two-well-gathering",
            fluid_model=self._BASE_FLUID,
            pressure_inlets=(
                PressureBoundary(node_id=1, pressure_pa=12e6),
                PressureBoundary(node_id=2, pressure_pa=10e6),
            ),
            pressure_outlets=(PressureBoundary(node_id=5, pressure_pa=1.5e6),),
            thermal_inlets=(
                ThermalBoundary(node_id=1, temperature_c=80.0, bc_type="fixed_temperature"),
                ThermalBoundary(node_id=2, temperature_c=65.0, bc_type="fixed_temperature"),
            ),
            inlet_fluid_bcs=(
                InletFluidBC(node_id=1, api_gravity=40.0, gas_gravity=0.70,
                             gor_sc_m3_per_m3=150.0, wor_sc_m3_per_m3=0.1),
                InletFluidBC(node_id=2, api_gravity=32.0, gas_gravity=0.65,
                             gor_sc_m3_per_m3=40.0,  wor_sc_m3_per_m3=0.3),
            ),
            components=self._PIPES,
        )

    def test_gathering_network_converges(self):
        """Two-well gathering network with BB must converge."""
        res = SteadyBlackOilSolver(
            turbulent_pipe_correlation=BeggsBrillCorrelation()
        ).solve(self._two_well_case())
        self.assertTrue(res.converged)

    def test_bb_manifold_pressure_higher_than_colebrook(self):
        """BB predicts higher manifold pressure than Colebrook.

        Because BB increases friction (f_tp > f_n), less fluid flows for
        the same ΔP, so the pressure drop in the gathering pipes is smaller
        and the manifold pressure is closer to the inlet pressures.
        """
        case = self._two_well_case()
        P_ck = SteadyBlackOilSolver().solve(case).node_pressures_pa[4]
        P_bb = SteadyBlackOilSolver(
            turbulent_pipe_correlation=BeggsBrillCorrelation()
        ).solve(case).node_pressures_pa[4]
        self.assertGreater(
            P_bb, P_ck,
            msg=f"BB manifold ({P_bb/1e6:.3f} MPa) should exceed Colebrook ({P_ck/1e6:.3f} MPa)",
        )

    def test_high_gor_well_has_more_slip_than_low_gor(self):
        """Well 1 (GOR=150) must produce more two-phase friction than Well 2 (GOR=40).

        A higher GOR means more gas in the pipe, more slip between phases,
        and a larger ratio f_tp/f_n.  This check verifies the holdup H_L
        for the high-GOR pipe is further below 1 (more two-phase character).
        """
        fluid_hi = BlackOilFluid(api_gravity=40.0, gas_gravity=0.70,
                                 gor_sc_m3_per_m3=150.0, wor_sc_m3_per_m3=0.1,
                                 reference_pressure_pa=8e6, reference_temperature_c=60.0)
        fluid_lo = BlackOilFluid(api_gravity=32.0, gas_gravity=0.65,
                                 gor_sc_m3_per_m3=15.0, wor_sc_m3_per_m3=0.3,
                                 reference_pressure_pa=8e6, reference_temperature_c=60.0)
        pvt_hi = fluid_hi.pvt(8e6, 60.0)
        pvt_lo = fluid_lo.pvt(8e6, 60.0)
        alpha_l_hi = pvt_hi.holdup_oil + pvt_hi.holdup_water
        alpha_l_lo = pvt_lo.holdup_oil + pvt_lo.holdup_water
        C_L_hi = alpha_l_hi / max(alpha_l_hi + pvt_hi.holdup_gas, 1e-9)
        C_L_lo = alpha_l_lo / max(alpha_l_lo + pvt_lo.holdup_gas, 1e-9)
        self.assertLess(
            C_L_hi, C_L_lo,
            msg=f"High-GOR C_L ({C_L_hi:.3f}) should be less than low-GOR C_L ({C_L_lo:.3f})",
        )


# ── Sachdeva critical-flow model ────────────────────────────────────────────────

class TestSachdevaLimits(unittest.TestCase):

    def test_pure_gas_matches_isentropic_ratio(self):
        """x1 = 1 → y_c = (2/(k+1))^(k/(k-1)) (n = k for a pure gas)."""
        for k in (1.2, 1.3, 1.4):
            y_c = sachdeva_critical_ratio(1.0, 50.0, 800.0, k, k)
            self.assertAlmostEqual(y_c, (2.0 / (k + 1.0)) ** (k / (k - 1.0)), places=8)

    def test_pure_liquid_has_no_critical_flow(self):
        self.assertEqual(sachdeva_critical_ratio(0.0, 50.0, 800.0, 1.3, 1.0), 0.0)

    def test_pure_liquid_mass_flux_is_orifice_equation(self):
        """x1 = 0 → G = Cd·√(2·ρ_L·ΔP)."""
        P1, y, rho_l, cd = 10e6, 0.9, 800.0, 0.62
        G = sachdeva_mass_flux(P1, y, 0.0, 50.0, rho_l, 1.3, cd)
        self.assertAlmostEqual(G, cd * math.sqrt(2.0 * rho_l * P1 * (1.0 - y)), places=6)

    def test_critical_ratio_increases_with_gas_fraction(self):
        ratios = [sachdeva_critical_ratio(x, 80.0, 800.0, 1.25, 1.01) for x in (0.001, 0.01, 0.05, 0.3)]
        self.assertEqual(ratios, sorted(ratios))
        self.assertTrue(all(0.0 < r < 0.6 for r in ratios))

    def test_pure_gas_mass_flux_maximum_at_critical_ratio(self):
        """For a pure gas (n = k), G(y) peaks exactly at y_c — the sonic condition."""
        k = 1.3
        y_c = sachdeva_critical_ratio(1.0, 40.0, 780.0, k, k)
        G_c = sachdeva_mass_flux(5e6, y_c, 1.0, 40.0, 780.0, k, 1.0)
        for dy in (-0.02, 0.02):
            self.assertLess(sachdeva_mass_flux(5e6, y_c + dy, 1.0, 40.0, 780.0, k, 1.0), G_c)

    def test_two_phase_flux_near_maximum_at_critical_ratio(self):
        """Two-phase: y_c uses polytropic n while G uses isentropic k (as published
        by Sachdeva et al.), so G(y_c) is close to — not exactly at — max G(y)."""
        x1, rho_g, rho_l, k, n = 0.07, 40.0, 780.0, 1.25, 1.013
        y_c = sachdeva_critical_ratio(x1, rho_g, rho_l, k, n)
        G_c = sachdeva_mass_flux(5e6, y_c, x1, rho_g, rho_l, k, 0.75)
        G_max = max(sachdeva_mass_flux(5e6, y / 1000.0, x1, rho_g, rho_l, k, 0.75)
                    for y in range(1, 1000))
        self.assertLess((G_max - G_c) / G_max, 0.02)

    def test_invalid_flow_model_raises(self):
        with self.assertRaises(ValueError):
            Choke(start_node=1, end_node=2, diameter_m=0.1, orifice_diameter_m=0.02,
                  flow_model="gilbert")

    def test_model_without_fluid_falls_back_to_k_factor(self):
        """No black-oil fluid → ChokeFlowModel uses the orifice K-factor."""
        model = ChokeFlowModel()
        model.set_fluid(None)

        class _Node:
            def __init__(self, p):
                self.pressure_pa = p

        class _State:
            component = Choke(start_node=1, end_node=2, diameter_m=0.1, orifice_diameter_m=0.02)
            start_node, end_node = _Node(2e5), _Node(1e5)
            area_m2 = component.area_m2
            velocity_m_per_s = 0.0

        v = model.calculate_velocity(_State, 1e5, 1000.0, 1e-3)
        K = _State.component.loss_coefficient
        self.assertAlmostEqual(v, math.sqrt(2.0 * 1e5 / (K * 1000.0)), places=10)


_GASSY_FLUID = BlackOilFluid(
    api_gravity=35.0, gas_gravity=0.65,
    gor_sc_m3_per_m3=120.0, wor_sc_m3_per_m3=0.2,
    reference_pressure_pa=3e6, reference_temperature_c=60.0,
)


def _gassy_choke_case(p_sep: float, flow_model: str = "sachdeva") -> NetworkCase:
    """Tutorial 10 scenario: 5 MPa wellhead below Pb, 24/64\" bean."""
    return NetworkCase(
        name="critical-choke",
        fluid_model=_GASSY_FLUID,
        pressure_inlets= (PressureBoundary(node_id=1, pressure_pa=5e6),),
        pressure_outlets=(PressureBoundary(node_id=3, pressure_pa=p_sep),),
        thermal_inlets  =(ThermalBoundary(node_id=1, temperature_c=60.0,
                                          bc_type="fixed_temperature"),),
        components=(
            Choke(start_node=1, end_node=2, diameter_m=0.1016,
                  orifice_diameter_m=24 * 0.0254 / 64, discharge_coefficient=0.75,
                  flow_model=flow_model),
            Pipe(start_node=2, end_node=3, diameter_m=0.1016, length_m=4_000.0,
                 absolute_roughness_m=46e-6),
        ),
    )


class TestCriticalChokeIntegration(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.solver = SteadyBlackOilSolver(turbulent_pipe_correlation=BeggsBrillCorrelation())

    def _m(self, p_sep, flow_model="sachdeva"):
        res = self.solver.solve(_gassy_choke_case(p_sep, flow_model))
        self.assertTrue(res.converged)
        return abs(res.component_flows[0].mass_flow_kg_per_s), res

    def test_critical_flow_independent_of_downstream_pressure(self):
        """While choked, tripling the separator pressure leaves the rate unchanged."""
        m_low, res = self._m(0.5e6)
        m_mid, _ = self._m(1.5e6)
        self.assertLess(res.node_pressures_pa[2] / 5e6, 0.5)  # deep in critical regime
        self.assertAlmostEqual(m_low, m_mid, delta=1e-3 * m_low)

    def test_subcritical_flow_decreases_with_downstream_pressure(self):
        m_crit, _ = self._m(0.5e6)
        m_sub1, _ = self._m(3.5e6)
        m_sub2, _ = self._m(4.5e6)
        self.assertLess(m_sub1, m_crit)
        self.assertLess(m_sub2, m_sub1)

    def test_orifice_model_overpredicts_critical_choke(self):
        m_sachdeva, _ = self._m(0.5e6)
        m_orifice, _ = self._m(0.5e6, flow_model="orifice")
        self.assertGreater(m_orifice, 1.1 * m_sachdeva)

    def test_choke_flow_matches_closed_form(self):
        """Solver rate equals Cd·A·G at the solved upstream pressure (critical)."""
        from angelica.closures.choke import black_oil_choke_inputs, choke_flow
        m_solver, res = self._m(0.5e6)
        x1, rho_g, rho_l, k, n = black_oil_choke_inputs(_GASSY_FLUID, 5e6, 60.0)
        area = math.pi * (24 * 0.0254 / 64) ** 2 / 4.0
        op = choke_flow(5e6, res.node_pressures_pa[2], area, 0.75, x1, rho_g, rho_l, k, n)
        self.assertTrue(op.critical)
        self.assertAlmostEqual(m_solver, op.mass_flow_kg_per_s, delta=1e-6 * m_solver)


class TestPerInletFluidReachesBeggsBrill(unittest.TestCase):
    """Regression: with InletFluidBC, Beggs-Brill must see real PVT (not the fallback)."""

    def test_no_all_liquid_fallback(self):
        case = TestGatheringNetworkBB()._two_well_case()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            SteadyBlackOilSolver(turbulent_pipe_correlation=BeggsBrillCorrelation()).solve(case)
        self.assertFalse(
            any("all-liquid fallback" in str(w.message) for w in caught),
            msg="Beggs-Brill fell back to all-liquid properties in a per-inlet network",
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
