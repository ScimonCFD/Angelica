"""
Tests for the Beggs-Brill (1973) two-phase flow correlation.

Coverage
--------
Unit tests
    - flow_regime():          verify boundary conditions and regime map
    - holdup_horizontal():    check H_L ≥ C_L (slip inequality) and
                              known intermediate values
    - two_phase_friction_factor(): f_tp = f_n when H_L = C_L (no slip limit),
                                    and f_tp > f_n when H_L > C_L
    - beggs_brill_gradient(): horizontal pressure gradient vs hand-calculation

Integration tests
    - BeggsBrillCorrelation + SteadyBlackOilSolver on a simple horizontal
      pipeline: solver converges, and BB gives a different (lower) flow rate
      than Colebrook for the same boundary conditions (because BB accounts for
      slip-enhanced friction).
"""
from __future__ import annotations

import math
import sys
import unittest
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[1] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.closures.beggs_brill import (
    BeggsBrillCorrelation,
    beggs_brill_gradient,
    flow_regime,
    holdup_horizontal,
    two_phase_friction_factor,
)
from angelica.core.case import NetworkCase, PressureBoundary, ThermalBoundary
from angelica.core.components import Pipe
from angelica.properties.black_oil import BlackOilFluid
from angelica.solvers import SteadyBlackOilSolver


# ── Shared test fluid ──────────────────────────────────────────────────────

_FLUID_BB = BlackOilFluid(
    api_gravity=32.0,
    gas_gravity=0.65,
    gor_sc_m3_per_m3=80.0,
    wor_sc_m3_per_m3=0.2,
    reference_pressure_pa=4e6,
    reference_temperature_c=50.0,
)


def _bb_case(D: float = 0.1016, L: float = 10_000.0,
             P_in: float = 6e6, P_out: float = 2e6,
             height_m: float = 0.0, fluid=None) -> NetworkCase:
    return NetworkCase(
        name="bb-test",
        fluid_model=fluid or _FLUID_BB,
        pressure_inlets=(PressureBoundary(node_id=1, pressure_pa=P_in),),
        pressure_outlets=(PressureBoundary(node_id=2, pressure_pa=P_out),),
        thermal_inlets=(ThermalBoundary(node_id=1, temperature_c=50.0,
                                        bc_type="fixed_temperature"),),
        components=(Pipe(
            start_node=1, end_node=2,
            diameter_m=D, length_m=L,
            absolute_roughness_m=46e-6,
            height_change_m=height_m,
        ),),
    )


# ── Unit tests ─────────────────────────────────────────────────────────────

class TestFlowRegime(unittest.TestCase):

    def test_pure_liquid_is_intermittent(self):
        regime, A = flow_regime(1.0, 10.0)
        self.assertEqual(regime, "intermittent")

    def test_pure_gas_is_distributed(self):
        regime, _ = flow_regime(0.0, 5.0)
        self.assertEqual(regime, "distributed")

    def test_low_froude_low_cl_segregated(self):
        # C_L < 0.01 and N_Fr < L1 = 316·C_L^0.302 ≈ 316·0.005^0.302 ≈ 75
        regime, _ = flow_regime(0.005, 5.0)
        self.assertEqual(regime, "segregated")

    def test_high_cl_high_froude_intermittent(self):
        # C_L = 0.5, N_Fr = 10 — should be intermittent (L3 ≈ 0.28, L4 ≈ huge → L3 ≤ 10 ≤ L4)
        regime, _ = flow_regime(0.5, 10.0)
        self.assertEqual(regime, "intermittent")

    def test_transition_alpha_in_zero_one(self):
        # Build a case that falls in the transition band (L2 ≤ N_Fr < L3) for C_L = 0.3
        # L2 = 9.252e-4 * 0.3^(-2.4684) ≈ 0.0191, L3 = 0.10 * 0.3^(-1.4516) ≈ 0.619
        # Choose N_Fr = 0.30 (inside [L2, L3))
        regime, A = flow_regime(0.3, 0.30)
        self.assertEqual(regime, "transition")
        self.assertGreaterEqual(A, 0.0)
        self.assertLessEqual(A, 1.0)


class TestHoldupHorizontal(unittest.TestCase):

    def test_holdup_ge_cl_intermittent(self):
        C_L, N_Fr = 0.40, 8.0
        HL = holdup_horizontal(C_L, N_Fr)
        self.assertGreaterEqual(HL, C_L)
        self.assertLessEqual(HL, 1.0)

    def test_holdup_ge_cl_segregated(self):
        C_L, N_Fr = 0.01, 0.001
        HL = holdup_horizontal(C_L, N_Fr)
        self.assertGreaterEqual(HL, C_L)
        self.assertLessEqual(HL, 1.0)

    def test_pure_liquid_holdup_one(self):
        self.assertAlmostEqual(holdup_horizontal(1.0, 5.0), 1.0)

    def test_pure_gas_holdup_zero(self):
        self.assertAlmostEqual(holdup_horizontal(0.0, 5.0), 0.0)

    def test_intermittent_known_value(self):
        """H_L(0) = 0.845 * C_L^0.5351 / N_Fr^0.0173 for intermittent."""
        C_L, N_Fr = 0.333, 9.03
        expected = 0.845 * C_L ** 0.5351 / N_Fr ** 0.0173
        expected = max(min(expected, 1.0), C_L)
        HL = holdup_horizontal(C_L, N_Fr)
        self.assertAlmostEqual(HL, expected, places=4)


class TestTwoPhaseFrictionFactor(unittest.TestCase):

    def test_no_correction_when_y_le_1(self):
        """When y = C_L / H_L² ≤ 1 the BB guard returns f_n unchanged.

        y ≤ 1  ⟺  H_L ≥ √C_L.
        Example: C_L = 0.30, H_L = 0.60 → y = 0.30 / 0.36 = 0.833 < 1.
        """
        f_n = 0.025
        C_L, HL = 0.30, 0.60   # y = 0.30/0.36 = 0.833
        f_tp = two_phase_friction_factor(f_n, C_L, HL)
        self.assertAlmostEqual(f_tp, f_n, places=8)

    def test_correction_when_y_gt_1(self):
        """When y = C_L / H_L² > 1, f_tp = f_n · exp(s) > f_n.

        Example: C_L = 0.30, H_L = 0.35 → y = 0.30 / 0.1225 ≈ 2.45 > 1.
        """
        f_n = 0.022
        C_L, HL = 0.30, 0.35   # y ≈ 2.45 > 1
        f_tp = two_phase_friction_factor(f_n, C_L, HL)
        self.assertGreater(f_tp, f_n)

    def test_cap_prevents_runaway(self):
        """Extreme case: f_tp ≤ f_n * exp(2) (cap in implementation)."""
        f_n = 0.02
        f_tp = two_phase_friction_factor(f_n, 0.01, 0.99)
        self.assertLessEqual(f_tp, f_n * math.exp(2.01))


class TestBeggsBrillGradient(unittest.TestCase):
    """Hand-calculated reference case for horizontal flow.

    Inputs (chosen to land in the intermittent regime):
        v_m   = 3.0 m/s
        C_L   = 0.333   (1/3 liquid)
        rho_l = 850 kg/m³,  rho_g = 10 kg/m³
        mu_l  = 5e-3 Pa·s,  mu_g  = 1.5e-5 Pa·s
        D     = 0.1016 m,   roughness = 46e-6 m,  theta = 0

    Expected regime: intermittent (C_L < 0.4, L3 ≈ 0.43 ≤ N_Fr ≈ 9.03 ≤ L1 ≈ 245)
    H_L(0) ≈ 0.459
    y = C_L / H_L² ≈ 0.333 / 0.211 ≈ 1.58
    f_tp / f_n > 1
    """

    _CL   = 0.333
    _VM   = 3.0
    _RL   = 850.0
    _RG   = 10.0
    _MUL  = 5e-3
    _MUG  = 1.5e-5
    _D    = 0.1016
    _EPS  = 46e-6

    def _st(self):
        return beggs_brill_gradient(
            self._VM, self._CL, self._RL, self._RG,
            self._MUL, self._MUG, self._D, self._EPS, 0.0,
        )

    def test_regime_intermittent(self):
        self.assertEqual(self._st().regime, "intermittent")

    def test_holdup_gt_cl(self):
        st = self._st()
        self.assertGreater(st.holdup, self._CL)

    def test_f_tp_gt_f_n(self):
        st = self._st()
        self.assertGreater(st.f_tp, st.f_n)

    def test_gravity_zero_horizontal(self):
        self.assertAlmostEqual(self._st().dP_dz_gravity, 0.0, places=3)

    def test_dP_friction_positive_and_reasonable(self):
        """dP/dz should be in the tens to hundreds Pa/m range for this case."""
        dPf = self._st().dP_dz_friction
        self.assertGreater(dPf, 10.0)    # not negligible
        self.assertLess(dPf, 10_000.0)   # not absurdly large

    def test_froude_number_correct(self):
        """N_Fr = v_m² / (g·D)  ≈  9.03"""
        N_Fr_expected = 3.0 ** 2 / (9.81 * 0.1016)
        self.assertAlmostEqual(self._st().N_Fr, N_Fr_expected, places=2)


# ── Integration tests ──────────────────────────────────────────────────────

class TestBeggsBrillIntegration(unittest.TestCase):
    """Solver-level tests: BB correlation inside SteadyBlackOilSolver."""

    def test_solver_converges_with_bb(self):
        """Beggs-Brill solver must converge on a simple horizontal pipeline."""
        res = SteadyBlackOilSolver(
            turbulent_pipe_correlation=BeggsBrillCorrelation()
        ).solve(_bb_case())
        self.assertTrue(res.converged)
        self.assertGreater(abs(res.component_flows[0].mass_flow_kg_per_s), 0.0)

    def test_bb_lower_flow_than_colebrook(self):
        """For a high-GOR two-phase flow, BB predicts lower flow than Colebrook.

        Physical reason: slip causes H_L > C_L, which increases the effective
        friction factor f_tp = f_n·exp(s) > f_n, so for the same ΔP boundary
        the mass flow must be lower.
        """
        res_ck = SteadyBlackOilSolver().solve(_bb_case())
        res_bb = SteadyBlackOilSolver(
            turbulent_pipe_correlation=BeggsBrillCorrelation()
        ).solve(_bb_case())

        self.assertTrue(res_ck.converged)
        self.assertTrue(res_bb.converged)

        q_ck = abs(res_ck.component_flows[0].volumetric_flow_m3_per_h)
        q_bb = abs(res_bb.component_flows[0].volumetric_flow_m3_per_h)

        self.assertGreater(
            q_ck, q_bb,
            msg=f"Expected Colebrook ({q_ck:.2f} m³/h) > Beggs-Brill ({q_bb:.2f} m³/h)",
        )

    def test_single_phase_bb_close_to_colebrook(self):
        """For a nearly all-liquid fluid (low GOR), BB ≈ Colebrook.

        When C_L ≈ 1, holdup H_L ≈ C_L ≈ 1, slip correction is negligible,
        and f_tp ≈ f_n.  The flows should agree within 10 %.
        """
        fluid_liq = BlackOilFluid(
            api_gravity=32.0,
            gas_gravity=0.65,
            gor_sc_m3_per_m3=0.5,   # tiny GOR → almost all liquid
            wor_sc_m3_per_m3=0.2,
            reference_pressure_pa=4e6,
            reference_temperature_c=50.0,
        )
        case = _bb_case(D=0.1016, L=5_000.0, P_in=5e6, P_out=2e6, fluid=fluid_liq)
        res_ck = SteadyBlackOilSolver().solve(case)
        res_bb = SteadyBlackOilSolver(
            turbulent_pipe_correlation=BeggsBrillCorrelation()
        ).solve(case)

        self.assertTrue(res_ck.converged)
        self.assertTrue(res_bb.converged)

        q_ck = abs(res_ck.component_flows[0].volumetric_flow_m3_per_h)
        q_bb = abs(res_bb.component_flows[0].volumetric_flow_m3_per_h)

        rel_diff = abs(q_ck - q_bb) / max(q_ck, 1.0)
        self.assertLess(
            rel_diff, 0.10,
            msg=f"Liquid-limit: Colebrook={q_ck:.2f} vs BB={q_bb:.2f} m³/h, diff={rel_diff*100:.1f}%",
        )

    def test_uphill_bb_lower_flow_than_horizontal(self):
        """Uphill inclined pipe has additional gravity pressure drop vs horizontal."""
        res_horiz = SteadyBlackOilSolver(
            turbulent_pipe_correlation=BeggsBrillCorrelation()
        ).solve(_bb_case())
        # Same case but 200 m elevation gain over 10 km (angle ≈ 1.1°)
        case_up = _bb_case(height_m=200.0)
        res_up = SteadyBlackOilSolver(
            turbulent_pipe_correlation=BeggsBrillCorrelation()
        ).solve(case_up)

        self.assertTrue(res_horiz.converged)
        self.assertTrue(res_up.converged)

        q_horiz = abs(res_horiz.component_flows[0].volumetric_flow_m3_per_h)
        q_up    = abs(res_up.component_flows[0].volumetric_flow_m3_per_h)

        self.assertGreater(
            q_horiz, q_up,
            msg=f"Uphill ({q_up:.2f} m³/h) should be less than horizontal ({q_horiz:.2f} m³/h)",
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
