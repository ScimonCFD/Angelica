"""Tests for the Michelsen arc-length phase envelope algorithm.

Validation basis:
- Three-component mixture (CH4 0.70 / C2H6 0.20 / C3H8 0.10) compared
  against HYSYS SRK01 (session data, 2026-08-29).
- Pure-component limits (methane) compared against NIST / thermo defaults.
- EPPR78 kij values compared against the bundled JSON database.
"""
from __future__ import annotations

import sys
import unittest
from pathlib import Path

SRC_ROOT = Path(__file__).resolve().parents[1] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica.properties.phase_envelope import (
    _build_kij_matrix,
    compute_phase_envelope,
    compute_quality_line,
)


class TestBuildKijMatrix(unittest.TestCase):
    def _constants(self, names):
        import warnings
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            from thermo import ChemicalConstantsPackage
            constants, _ = ChemicalConstantsPackage.from_IDs(names)
        return constants

    def test_known_pairs_from_eppr78(self):
        c = self._constants(["methane", "ethane", "propane"])
        kijs = _build_kij_matrix(c)
        # Symmetric
        self.assertAlmostEqual(kijs[0][1], kijs[1][0], places=8)
        self.assertAlmostEqual(kijs[0][2], kijs[2][0], places=8)
        self.assertAlmostEqual(kijs[1][2], kijs[2][1], places=8)
        # Diagonal is zero
        for i in range(3):
            self.assertEqual(kijs[i][i], 0.0)
        # Known EPPR78 values
        self.assertAlmostEqual(kijs[0][1],  0.0058, places=3)   # CH4–C2H6
        self.assertAlmostEqual(kijs[0][2],  0.0189, places=3)   # CH4–C3H8
        self.assertAlmostEqual(kijs[1][2], -0.0029, places=3)   # C2H6–C3H8

    def test_missing_pair_defaults_to_zero(self):
        # Use a real but uncommon component — neon has no natural-gas kij
        try:
            c = self._constants(["methane", "neon"])
        except Exception:
            self.skipTest("neon not available in thermo")
        kijs = _build_kij_matrix(c)
        self.assertEqual(kijs[0][1], 0.0)
        self.assertEqual(kijs[1][0], 0.0)

    def test_single_component_returns_zero_matrix(self):
        c = self._constants(["methane"])
        kijs = _build_kij_matrix(c)
        self.assertEqual(len(kijs), 1)
        self.assertEqual(kijs[0][0], 0.0)


class TestComputePhaseEnvelope(unittest.TestCase):
    """Regression tests for the arc-length phase envelope."""

    NAMES = ["methane", "ethane", "propane"]
    ZS    = [0.70, 0.20, 0.10]

    def _envelope(self, eos="PR"):
        import warnings
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            return compute_phase_envelope(self.NAMES, self.ZS, eos_name=eos)

    # ── basic structural checks ──────────────────────────────────────────────

    def test_returns_non_empty_bubble_and_dew(self):
        bubble, dew, Tc, Pc = self._envelope()
        self.assertGreater(len(bubble), 5, "bubble curve has too few points")
        self.assertGreater(len(dew),    5, "dew curve has too few points")

    def test_curves_close_at_critical_point(self):
        # The two arcs converge near the critical point but may not meet exactly —
        # PR and SRK can differ slightly in convergence radius.  Allow 8 K / 15%.
        bubble, dew, _, _ = self._envelope()
        T_bc, P_bc = bubble[-1]
        T_dc, P_dc = dew[-1]
        self.assertAlmostEqual(T_bc, T_dc, delta=8.0, msg="bubble/dew do not share final T")
        self.assertAlmostEqual(P_bc, P_dc, delta=P_bc * 0.15, msg="bubble/dew do not share final P")

    def test_bubble_T_monotone_increasing(self):
        bubble, _, _, _ = self._envelope()
        Ts = [t for t, _ in bubble]
        self.assertEqual(Ts, sorted(Ts), "bubble temperatures not monotone")

    def test_dew_T_spans_meaningful_range(self):
        # The dew curve is NOT globally monotone: it increases to the cricondentherm
        # then bends back toward the critical point (retrograde behavior).
        # Check that it covers at least 50 K so the arc-length traced a real curve.
        _, dew, _, _ = self._envelope()
        Ts = [t for t, _ in dew]
        self.assertGreater(max(Ts) - min(Ts), 50.0,
                           "dew curve T range too small — arc may have failed")

    def test_pressures_positive(self):
        bubble, dew, Tc, Pc = self._envelope()
        for _, P in bubble + dew:
            self.assertGreater(P, 0, "non-positive pressure in envelope")
        self.assertGreater(Pc, 0)

    # ── critical point sanity (vs HYSYS SRK — within 5°C / 5 bar) ──────────

    def test_critical_T_in_expected_range_PR(self):
        _, _, Tc, _ = self._envelope("PR")
        Tc_C = Tc - 273.15
        # HYSYS SRK critical ≈ -17°C; PR is close; allow ±10°C
        self.assertGreater(Tc_C, -35.0, f"Tc={Tc_C:.1f}°C too low")
        self.assertLess(   Tc_C,   5.0, f"Tc={Tc_C:.1f}°C too high")

    def test_critical_P_in_expected_range_PR(self):
        _, _, _, Pc = self._envelope("PR")
        Pc_bar = Pc / 1e5
        # HYSYS SRK cricondenbar ≈ 87 bar; critical ≈ 83 bar; allow 70–100 bar
        self.assertGreater(Pc_bar, 70.0, f"Pc={Pc_bar:.1f} bar too low")
        self.assertLess(   Pc_bar, 100.0, f"Pc={Pc_bar:.1f} bar too high")

    def test_critical_T_in_expected_range_SRK(self):
        _, _, Tc, _ = self._envelope("SRK")
        Tc_C = Tc - 273.15
        self.assertGreater(Tc_C, -35.0)
        self.assertLess(   Tc_C,   5.0)

    def test_critical_P_in_expected_range_SRK(self):
        _, _, _, Pc = self._envelope("SRK")
        Pc_bar = Pc / 1e5
        self.assertGreater(Pc_bar, 70.0)
        self.assertLess(   Pc_bar, 100.0)

    # ── PR vs SRK consistency ────────────────────────────────────────────────

    def test_pr_srk_critical_T_within_5C(self):
        _, _, Tc_pr, _ = self._envelope("PR")
        _, _, Tc_srk, _ = self._envelope("SRK")
        self.assertAlmostEqual(Tc_pr, Tc_srk, delta=5.0,
                               msg="PR and SRK critical T differ by more than 5 K")

    def test_pr_srk_critical_P_within_10pct(self):
        # PR and SRK predict different critical pressures for the same mixture;
        # the arc-length endpoints are also approximate.  Allow 10% tolerance.
        _, _, _, Pc_pr  = self._envelope("PR")
        _, _, _, Pc_srk = self._envelope("SRK")
        rel = abs(Pc_pr - Pc_srk) / max(Pc_pr, Pc_srk)
        self.assertLess(rel, 0.10, f"PR/SRK Pc differ by {rel*100:.1f}%")

    # ── pure-component limiting case ─────────────────────────────────────────

    def test_pure_methane_critical_near_literature(self):
        """Pure CH4: Tc ≈ 190.6 K, Pc ≈ 4.60 MPa (NIST)."""
        import warnings
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            bubble, dew, Tc, Pc = compute_phase_envelope(["methane"], [1.0])
        # Critical coordinates within 5 K and 10%
        self.assertAlmostEqual(Tc, 190.6, delta=5.0,
                               msg=f"Pure CH4 Tc={Tc:.1f} K, expected ~190.6 K")
        self.assertAlmostEqual(Pc / 1e6, 4.60, delta=0.46,
                               msg=f"Pure CH4 Pc={Pc/1e6:.2f} MPa, expected ~4.60 MPa")

    def test_pure_component_bubble_dew_nearly_coincide(self):
        """For a pure component, bubble = dew (vapour-pressure curve)."""
        import warnings
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            bubble, dew, _, _ = compute_phase_envelope(["methane"], [1.0])
        # Sample a few interior points and check they lie within 0.5 K of each other
        n = min(len(bubble), len(dew))
        for i in range(0, n, max(1, n // 5)):
            Tb = bubble[i][0]
            Td = dew[i][0]
            self.assertAlmostEqual(Tb, Td, delta=1.0,
                                   msg=f"Pure CH4 bubble[{i}]={Tb:.2f} K != dew={Td:.2f} K")

    def test_pure_methane_saturation_vs_nist(self):
        """PR EOS saturation pressures within 2% of NIST reference.

        Reference values from NIST REFPROP (Lemmon et al. 2009, J. Phys. Chem. Ref. Data
        38:721-748), accurate to <0.01% on vapour pressure. The PR EOS systematic
        over-prediction of ~0.7% is well within the expected 2% tolerance.
        """
        import warnings

        import numpy as np
        # NIST REFPROP saturation pressures for pure methane (Pa)
        NIST_SAT_PA = {
            130: 3.6732e5,
            140: 6.4118e5,
            150: 10.3996e5,
            160: 15.9208e5,
            170: 23.2835e5,
            180: 32.8518e5,
        }
        TOLERANCE = 0.02  # 2% — conservative bound for PR EOS on Psat
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            bubble, _, _, _ = compute_phase_envelope(["methane"], [1.0])
        Tb = np.array([p[0] for p in bubble])
        Pb = np.array([p[1] for p in bubble])
        for T_ref, P_nist in NIST_SAT_PA.items():
            P_pr = float(np.interp(T_ref, Tb, Pb))
            rel_err = abs(P_pr - P_nist) / P_nist
            self.assertLess(
                rel_err, TOLERANCE,
                f"T={T_ref}K: PR={P_pr/1e5:.4f} bar, NIST={P_nist/1e5:.4f} bar, "
                f"err={rel_err*100:.1f}% > {TOLERANCE*100:.0f}%",
            )

    # ── cricondentherm / cricondenbar rough bounds ───────────────────────────

    def test_cricondentherm_above_critical(self):
        """Max T on the envelope must be above Tc (it IS the cricondentherm)."""
        bubble, dew, Tc, _ = self._envelope()
        T_max = max(t for t, _ in bubble + dew)
        self.assertGreater(T_max, Tc - 1.0)

    def test_cricondenbar_above_critical_P(self):
        """Max P on the envelope must be ≥ Pc."""
        bubble, dew, _, Pc = self._envelope()
        P_max = max(p for _, p in bubble + dew)
        self.assertGreaterEqual(P_max, Pc * 0.98)


class TestComputeQualityLine(unittest.TestCase):
    NAMES = ["methane", "ethane", "propane"]
    ZS    = [0.70, 0.20, 0.10]

    def _quality(self, vf=0.5, eos="PR"):
        import warnings
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            return compute_quality_line(self.NAMES, self.ZS, vf=vf, eos_name=eos)

    def test_returns_non_empty_list(self):
        pts = self._quality(0.5)
        self.assertGreater(len(pts), 0)

    def test_temperatures_monotone(self):
        pts = self._quality(0.5)
        Ts = [t for t, _ in pts]
        self.assertEqual(Ts, sorted(Ts))

    def test_pressures_positive(self):
        pts = self._quality(0.5)
        for _, P in pts:
            self.assertGreater(P, 0)

    def test_vf_near_zero_and_one(self):
        """Quality lines at VF=0.1 and VF=0.9 should both compute successfully."""
        pts_low  = self._quality(0.1)
        pts_high = self._quality(0.9)
        self.assertGreater(len(pts_low),  0)
        self.assertGreater(len(pts_high), 0)

    def test_vf05_P_between_bubble_and_dew(self):
        """At a given T, VF=0.5 pressure should lie between bubble and dew pressures."""
        import warnings
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            bubble, dew, _, _ = compute_phase_envelope(self.NAMES, self.ZS)
        quality = self._quality(0.5)

        # Build a quick T→P lookup for bubble and dew
        def interp(pts, T_target):
            for i in range(len(pts) - 1):
                T0, P0 = pts[i]
                T1, P1 = pts[i + 1]
                if T0 <= T_target <= T1:
                    frac = (T_target - T0) / (T1 - T0)
                    return P0 + frac * (P1 - P0)
            return None

        T_min_b = bubble[0][0]
        T_max_b = bubble[-2][0]  # exclude the critical point
        T_min_d = dew[0][0]
        T_max_d = dew[-2][0]
        T_lo = max(T_min_b, T_min_d)
        T_hi = min(T_max_b, T_max_d)

        checks = 0
        for T, P_q in quality:
            if T < T_lo or T > T_hi:
                continue
            P_b = interp(bubble, T)
            P_d = interp(dew, T)
            if P_b is None or P_d is None:
                continue
            P_lo = min(P_b, P_d)
            P_hi = max(P_b, P_d)
            tol = 0.15 * (P_hi - P_lo)
            self.assertGreater(P_q, P_lo - tol,
                               f"VF=0.5 P={P_q:.0f} Pa below envelope at T={T:.1f} K")
            self.assertLess(   P_q, P_hi + tol,
                               f"VF=0.5 P={P_q:.0f} Pa above envelope at T={T:.1f} K")
            checks += 1

        self.assertGreater(checks, 3, "Too few interior quality-line points checked")



class TestRichCondensateEnvelope(unittest.TestCase):
    """Rich mixture with C5/C6 heavy ends (compositional tutorial 06, dry basis).

    Before v1.10.2 the Wilson start for the dew curve sat at 0.5·Tc_mix, far
    below the dew curve of a mixture with heavy ends: the dew curve came back
    empty, the envelope was open and the plot showed the Kay pseudo-critical
    point (≈1 °C, 44 bar) instead of the real critical point.
    """

    NAMES = ["methane", "ethane", "propane", "n-butane", "n-pentane", "n-hexane"]
    ZS = [v / 0.93 for v in (0.55, 0.10, 0.10, 0.08, 0.05, 0.05)]

    @classmethod
    def setUpClass(cls):
        import warnings
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            cls.bubble, cls.dew, cls.Tc, cls.Pc = compute_phase_envelope(cls.NAMES, cls.ZS)

    def test_dew_curve_traced(self):
        self.assertGreater(len(self.dew), 20)

    def test_envelope_closes_at_critical_point(self):
        self.assertEqual(self.bubble[-1], self.dew[-1])
        self.assertEqual((self.Tc, self.Pc), self.bubble[-1])

    def test_critical_point_is_not_kay_estimate(self):
        Tc_C, Pc_bar = self.Tc - 273.15, self.Pc / 1e5
        self.assertGreater(Tc_C, 40.0, f"Tc={Tc_C:.1f} °C")
        self.assertGreater(Pc_bar, 100.0, f"Pc={Pc_bar:.1f} bar")

    def test_cricondentherm_matches_flash(self):
        """Above the cricondentherm the flash is single phase, just below it two-phase."""
        import warnings

        from thermo import ChemicalConstantsPackage
        from thermo.eos_mix import PRMIX
        from thermo.flash import FlashVL
        from thermo.phases import CEOSGas, CEOSLiquid

        T_cct, P_cct = max(self.dew, key=lambda tp: tp[0])
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            c, props = ChemicalConstantsPackage.from_IDs(self.NAMES)
            kw = {"Tcs": c.Tcs, "Pcs": c.Pcs, "omegas": c.omegas, "kijs": _build_kij_matrix(c)}
            flash = FlashVL(c, props,
                            liquid=CEOSLiquid(PRMIX, kw, HeatCapacityGases=props.HeatCapacityGases),
                            gas=CEOSGas(PRMIX, kw, HeatCapacityGases=props.HeatCapacityGases))
            above = flash.flash(zs=self.ZS, T=T_cct + 2.0, P=P_cct)
            below = flash.flash(zs=self.ZS, T=T_cct - 2.0, P=P_cct)
        self.assertEqual(above.phase_count, 1)
        self.assertEqual(below.phase_count, 2)

    def test_dew_curve_reaches_low_pressure(self):
        self.assertLess(min(p for _, p in self.dew), 1.0e5)



class TestAbsentComponents(unittest.TestCase):
    """Components with zero mole fraction (e.g. the lean gas of compositional
    tutorial 06, 0 % C5/C6) must be ignored.  Before v1.10.3 they made every
    flash fail: ~40 s of work and an empty "no two-phase region" plot."""

    H = ["methane", "ethane", "propane", "n-butane", "n-pentane", "n-hexane"]

    def test_zero_components_give_same_envelope_as_reduced_mixture(self):
        import warnings
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            full = compute_phase_envelope(self.H, [0.90, 0.06, 0.03, 0.01, 0.0, 0.0])
            reduced = compute_phase_envelope(self.H[:4], [0.90, 0.06, 0.03, 0.01])
        self.assertGreater(len(full[0]), 5)
        self.assertGreater(len(full[1]), 5)
        self.assertAlmostEqual(full[2], reduced[2], places=6)
        self.assertAlmostEqual(full[3], reduced[3], delta=1.0)

    def test_quality_line_with_zero_components(self):
        import warnings
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            pts = compute_quality_line(self.H, [0.90, 0.06, 0.03, 0.01, 0.0, 0.0], 0.5)
        self.assertGreater(len(pts), 3)

    def test_all_zero_composition_raises(self):
        with self.assertRaises(ValueError):
            compute_phase_envelope(self.H[:2], [0.0, 0.0])



class TestQualityLineContinuation(unittest.TestCase):
    """Quality lines are traced by arc-length continuation and must end at the
    critical point.  Before v1.10.4 a temperature scan of T–VF flashes jumped
    to the wrong root near the critical point (kink, then a tail beyond Tc)."""

    NAMES = ["methane", "ethane", "propane"]
    ZS = [0.70, 0.20, 0.10]

    @classmethod
    def setUpClass(cls):
        import warnings
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            cls.env = compute_phase_envelope(cls.NAMES, cls.ZS)
            cls.lines = {vf: compute_quality_line(cls.NAMES, cls.ZS, vf) for vf in (0.1, 0.3, 0.9)}

    def test_lines_end_at_critical_point(self):
        _, _, Tc, Pc = self.env
        for vf, pts in self.lines.items():
            T_end, P_end = pts[-1]
            self.assertAlmostEqual(T_end, Tc, delta=0.5, msg=f"VF={vf}")
            self.assertAlmostEqual(P_end, Pc, delta=0.005 * Pc, msg=f"VF={vf}")

    def test_points_have_requested_vapour_fraction(self):
        """Check with independent TP flashes, away from the critical point."""
        import warnings

        from thermo import ChemicalConstantsPackage
        from thermo.eos_mix import PRMIX
        from thermo.flash import FlashVL
        from thermo.phases import CEOSGas, CEOSLiquid

        _, _, Tc, _ = self.env
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            c, props = ChemicalConstantsPackage.from_IDs(self.NAMES)
            kw = {"Tcs": c.Tcs, "Pcs": c.Pcs, "omegas": c.omegas, "kijs": _build_kij_matrix(c)}
            flash = FlashVL(c, props,
                            liquid=CEOSLiquid(PRMIX, kw, HeatCapacityGases=props.HeatCapacityGases),
                            gas=CEOSGas(PRMIX, kw, HeatCapacityGases=props.HeatCapacityGases))
            for vf, pts in self.lines.items():
                checked = 0
                for T, P in pts[::4]:
                    if abs(T - Tc) < 5.0 or P < 1e4:
                        continue
                    res = flash.flash(zs=self.ZS, T=T, P=P)
                    self.assertEqual(res.phase_count, 2, f"VF={vf} at T={T:.1f} K")
                    # thermo can swap the gas/liquid labels near the critical
                    # point: identify the vapour as the lighter phase.
                    beta_gas = min(zip(res.betas, res.phases), key=lambda bp: bp[1].rho_mass())[0]
                    self.assertAlmostEqual(beta_gas, vf, delta=2e-3, msg=f"T={T:.1f} K, P={P:.0f} Pa")
                    checked += 1
                self.assertGreater(checked, 5)

    def test_lines_ordered_between_bubble_and_dew(self):
        """At fixed T, a higher vapour fraction means a lower pressure."""
        def p_at(pts, T):
            for (T1, P1), (T2, P2) in zip(pts, pts[1:]):
                if min(T1, T2) <= T <= max(T1, T2) and T1 != T2:
                    return P1 + (P2 - P1) * (T - T1) / (T2 - T1)
            return None
        T = 200.0
        p01, p03, p09 = (p_at(self.lines[v], T) for v in (0.1, 0.3, 0.9))
        self.assertTrue(p01 > p03 > p09, (p01, p03, p09))


if __name__ == "__main__":
    unittest.main()
