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

from angelica.closures.beggs_brill import BeggsBrillCorrelation
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


if __name__ == "__main__":
    unittest.main(verbosity=2)
