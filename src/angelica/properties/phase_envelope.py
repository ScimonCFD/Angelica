"""Phase envelope (bubble/dew point curves) computation for multi-component mixtures.

Uses thermo.flash.FlashVL with Peng-Robinson EOS.  FlashVL properly fails
(raises or returns None) when bubble/dew conditions don't exist above the
mixture critical temperature, unlike the older thermo.Mixture(VF=) API which
returns spurious supercritical solutions.
"""

from __future__ import annotations


def compute_phase_envelope(
    component_names: tuple[str, ...] | list[str],
    zs: tuple[float, ...] | list[float],
    n_T: int = 20,
) -> tuple[list[tuple[float, float]], list[tuple[float, float]], float, float]:
    """Compute bubble and dew point curves for a multi-component mixture.

    Returns:
        (bubble_pts, dew_pts, Tc_K, Pc_Pa) ordered low→high temperature.
        Both curves share the same last point (the critical point) so the
        envelope is closed.  Tc_K / Pc_Pa are mole-fraction-weighted
        pseudo-critical coordinates used only for the plot marker.
    """
    import numpy as np
    from thermo import ChemicalConstantsPackage, PropertyCorrelationsPackage
    from thermo.flash import FlashVL
    from thermo.phases import CEOSGas, CEOSLiquid
    from thermo.eos_mix import PRMIX

    names = list(component_names)
    fracs = list(zs)

    constants, props = ChemicalConstantsPackage.from_IDs(names)

    eos_kw = dict(Tcs=constants.Tcs, Pcs=constants.Pcs, omegas=constants.omegas)
    gas_phase = CEOSGas(PRMIX, eos_kw, HeatCapacityGases=props.HeatCapacityGases)
    liq_phase = CEOSLiquid(PRMIX, eos_kw, HeatCapacityGases=props.HeatCapacityGases)
    flash_obj = FlashVL(constants, props, liquid=liq_phase, gas=gas_phase)

    Tcs: list[float] = list(constants.Tcs)
    Pcs: list[float] = list(constants.Pcs)
    Tbs: list[float] = list(constants.Tbs)

    Tc = sum(z * tc for z, tc in zip(fracs, Tcs))
    Pc = sum(z * pc for z, pc in zip(fracs, Pcs))

    T_lo = max(min(Tbs) * 0.65, 80.0)
    # The true mixture critical temperature (from the EOS) can exceed the
    # mole-fraction-weighted pseudo-Tc by 10–20 % for lean natural gas.
    # Scan to 1.35×Tc_pseudo to ensure we always cover the critical region.
    T_hi_scan = Tc * 1.35
    if T_lo >= T_hi_scan:
        T_lo = T_hi_scan * 0.50

    # ── Flash helpers ──────────────────────────────────────────────────────────
    # FlashVL correctly raises / returns a degenerate state above the mixture
    # critical temperature, so bubble and dew calculations naturally stop there.

    def _bubble(T: float) -> float | None:
        try:
            res = flash_obj.flash(zs=fracs, T=float(T), VF=0)
            P = res.P
            return float(P) if P is not None and P > 0 else None
        except Exception:
            return None

    def _dew(T: float) -> float | None:
        try:
            res = flash_obj.flash(zs=fracs, T=float(T), VF=1)
            P = res.P
            return float(P) if P is not None and P > 0 else None
        except Exception:
            return None

    # ── Temperature scan ───────────────────────────────────────────────────────
    T_vals = np.linspace(T_lo, T_hi_scan, n_T + 10)

    bubble_pts: list[tuple[float, float]] = []
    dew_pts: list[tuple[float, float]] = []
    envelope_closed = False
    last_both: tuple[float, float, float] | None = None  # (T, P_bub, P_dew)

    for T in T_vals:
        T = float(T)
        P_b = _bubble(T)
        P_d = _dew(T)

        if P_b is not None and P_d is not None:
            frac_diff = abs(P_b - P_d) / max(P_b, P_d)
            if frac_diff < 0.05:
                P_crit = (P_b + P_d) / 2.0
                bubble_pts.append((T, P_crit))
                dew_pts.append((T, P_crit))
                envelope_closed = True
                break
            bubble_pts.append((T, P_b))
            dew_pts.append((T, P_d))
            last_both = (T, P_b, P_d)
        elif P_b is not None:
            bubble_pts.append((T, P_b))
        elif P_d is not None:
            dew_pts.append((T, P_d))
        elif bubble_pts or dew_pts:
            break  # both failed — past the two-phase region

    # ── Closing bisection ─────────────────────────────────────────────────────
    # Bisect between the last temperature where both curves coexisted and the
    # first temperature where one of them failed, to find the precise critical
    # point (where P_bub ≈ P_dew).
    if not envelope_closed and last_both is not None:
        T_cl_lo, _, _ = last_both
        step = (T_hi_scan - T_lo) / (n_T + 10 - 1)
        T_cl_hi = T_cl_lo + step * 4

        best_T: float | None = None
        best_P: float | None = None
        best_frac = float("inf")

        for _ in range(20):
            if T_cl_hi - T_cl_lo < 0.05:
                break
            T_mid = (T_cl_lo + T_cl_hi) / 2.0
            P_b = _bubble(T_mid)
            P_d = _dew(T_mid)

            if P_b is None or P_d is None:
                T_cl_hi = T_mid
                continue

            frac_diff = abs(P_b - P_d) / max(P_b, P_d)
            if frac_diff < best_frac:
                best_frac = frac_diff
                best_T = T_mid
                best_P = (P_b + P_d) / 2.0

            if frac_diff < 0.05:
                break
            T_cl_lo = T_mid

        if best_T is not None and best_P is not None:
            bubble_pts = [(t, p) for t, p in bubble_pts if t < best_T]
            dew_pts = [(t, p) for t, p in dew_pts if t < best_T]
            bubble_pts.append((best_T, best_P))
            dew_pts.append((best_T, best_P))
        else:
            # Force-close at the last T where both curves coexisted.
            T_fc, P_fb, P_fd = last_both
            P_avg = (P_fb + P_fd) / 2.0
            bubble_pts = [(t, p) for t, p in bubble_pts if t <= T_fc]
            dew_pts = [(t, p) for t, p in dew_pts if t <= T_fc]
            if not bubble_pts or bubble_pts[-1][0] < T_fc:
                bubble_pts.append((T_fc, P_avg))
            if not dew_pts or dew_pts[-1][0] < T_fc:
                dew_pts.append((T_fc, P_avg))

    return bubble_pts, dew_pts, Tc, Pc
