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
    from thermo import ChemicalConstantsPackage
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

    # ── Temperature scan with hot_start ───────────────────────────────────────
    # Using hot_start passes the previous flash result as the initial K-factor
    # guess for the next temperature.  This dramatically stabilises convergence
    # near the critical point (the same technique used by thermo's plot_TP with
    # hot=True) and avoids the spurious solutions that occur with cold starts.
    T_vals = np.linspace(T_lo, T_hi_scan, n_T + 10)

    bubble_pts: list[tuple[float, float]] = []
    dew_pts: list[tuple[float, float]] = []
    envelope_closed = False
    last_both: tuple[float, float, float] | None = None  # (T, P_bub, P_dew)

    state_bub = None  # hot_start state for bubble branch
    state_dew = None  # hot_start state for dew branch

    for T in T_vals:
        T = float(T)

        P_b: float | None = None
        try:
            res = flash_obj.flash(zs=fracs, T=T, VF=0, hot_start=state_bub)
            if res.P is not None and res.P > 0:
                P_b = float(res.P)
                state_bub = res
            else:
                state_bub = None
        except Exception:
            state_bub = None

        P_d: float | None = None
        try:
            res = flash_obj.flash(zs=fracs, T=T, VF=1, hot_start=state_dew)
            if res.P is not None and res.P > 0:
                P_d = float(res.P)
                state_dew = res
            else:
                state_dew = None
        except Exception:
            state_dew = None

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

    # ── Fine-scan closing ─────────────────────────────────────────────────────
    # After the coarse scan breaks, step forward at 0.5 K intervals to fill in
    # the closing region.  Cold-start flashes are used here because near the
    # critical point they reliably find the convergent (K_i → 1) solution that
    # marks the closure, while hot-start can get trapped on a spurious branch.
    if not envelope_closed and last_both is not None:
        T_fc, _, _ = last_both
        fine_step = 0.1
        T_fine = T_fc + fine_step

        while T_fine <= T_hi_scan + fine_step:
            T_f = float(T_fine)
            P_b: float | None = None
            P_d: float | None = None
            try:
                rb = flash_obj.flash(zs=fracs, T=T_f, VF=0)
                if rb.P is not None and rb.P > 0:
                    P_b = float(rb.P)
            except Exception:
                pass
            try:
                rd = flash_obj.flash(zs=fracs, T=T_f, VF=1)
                if rd.P is not None and rd.P > 0:
                    P_d = float(rd.P)
            except Exception:
                pass

            if P_b is not None and P_d is not None:
                frac_diff = abs(P_b - P_d) / max(P_b, P_d)
                if frac_diff < 0.05:
                    P_crit = (P_b + P_d) / 2.0
                    bubble_pts.append((T_f, P_crit))
                    dew_pts.append((T_f, P_crit))
                    envelope_closed = True
                    break
                bubble_pts.append((T_f, P_b))
                dew_pts.append((T_f, P_d))
                last_both = (T_f, P_b, P_d)
            elif P_b is not None:
                bubble_pts.append((T_f, P_b))
            elif P_d is not None:
                dew_pts.append((T_f, P_d))
            else:
                break

            T_fine += fine_step

        if not envelope_closed and last_both is not None:
            T_fc2, P_fb, P_fd = last_both
            P_avg = (P_fb + P_fd) / 2.0
            bubble_pts = [(t, p) for t, p in bubble_pts if t <= T_fc2]
            dew_pts    = [(t, p) for t, p in dew_pts    if t <= T_fc2]
            if not bubble_pts or bubble_pts[-1][0] < T_fc2:
                bubble_pts.append((T_fc2, P_avg))
            if not dew_pts or dew_pts[-1][0] < T_fc2:
                dew_pts.append((T_fc2, P_avg))

    return bubble_pts, dew_pts, Tc, Pc
