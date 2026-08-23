"""Phase envelope (bubble/dew point curves) computation for multi-component mixtures.

Uses thermo.flash.FlashVL with Peng-Robinson EOS.  FlashVL properly fails
(raises or returns None) when bubble/dew conditions don't exist above the
mixture critical temperature, unlike the older thermo.Mixture(VF=) API which
returns spurious supercritical solutions.
"""

from __future__ import annotations


def _fill_cliff(
    bubble_pts: list[tuple[float, float]],
    T_close: float,
    P_close: float,
    n_fill: int = 8,
) -> None:
    """Insert square-root interpolated points when the closing step is a cliff.

    Applies P(T) = P_c + (P_last - P_c) * [(T_c - T)/(T_c - T_last)]^0.5,
    which matches the mean-field critical exponent of cubic EOS.  Only runs
    when the last bubble pressure exceeds the closing pressure by >30 %.
    """
    if not bubble_pts:
        return
    T_lb, P_lb = bubble_pts[-1]
    if P_lb <= P_close * 1.3 or T_close <= T_lb:
        return
    dT = T_close - T_lb
    for i in range(1, n_fill):
        frac = i / n_fill
        T_int = T_lb + dT * frac
        rem = 1.0 - frac
        P_int = P_close + (P_lb - P_close) * (rem ** 0.5)
        bubble_pts.append((T_int, P_int))


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

    # ── Temperature scan with hot_start ───────────────────────────────────────
    T_vals = np.linspace(T_lo, T_hi_scan, n_T + 10)

    bubble_pts: list[tuple[float, float]] = []
    dew_pts: list[tuple[float, float]] = []
    envelope_closed = False
    last_both: tuple[float, float, float] | None = None  # (T, P_bub_Pa, P_dew_Pa)

    state_bub = None  # hot_start state for bubble branch
    state_dew = None  # hot_start state for dew branch
    last_bub_state = None  # saved at the last T where both curves coexist
    last_dew_state = None

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
                _fill_cliff(bubble_pts, T, P_crit)
                bubble_pts.append((T, P_crit))
                dew_pts.append((T, P_crit))
                envelope_closed = True
                break
            bubble_pts.append((T, P_b))
            dew_pts.append((T, P_d))
            last_both = (T, P_b, P_d)
            last_bub_state = state_bub
            last_dew_state = state_dew
        elif P_b is not None:
            bubble_pts.append((T, P_b))
        elif P_d is not None:
            dew_pts.append((T, P_d))
        elif bubble_pts or dew_pts:
            break  # both failed — past the two-phase region

    # ── Fine-scan closing with hot_start ──────────────────────────────────────
    # Step forward at 0.5 K from the last coarse-scan point using hot_start.
    # hot_start keeps the bubble branch on the physical retrograde solution
    # (smoothly descending from the cricondenbar toward the critical point)
    # rather than jumping to the K_i→1 trivial solution that cold-start finds.
    # A spike filter rejects hot_start bubble results that rise more than 5 %
    # above the last accepted value, preventing isolated numerical instabilities.
    if not envelope_closed and last_both is not None:
        T_fc, _, _ = last_both
        fine_step = 0.5
        T_fine = T_fc + fine_step

        fine_bub_state = last_bub_state
        fine_dew_state = last_dew_state

        best_frac = float("inf")
        best_close_T: float | None = None
        best_close_P: float | None = None

        while T_fine <= T_hi_scan + fine_step:
            T_f = float(T_fine)

            # ── bubble ─────────────────────────────────────────────────────
            P_b: float | None = None
            prev_bub_state = fine_bub_state
            try:
                rb = flash_obj.flash(zs=fracs, T=T_f, VF=0, hot_start=fine_bub_state)
                if rb.P is not None and rb.P > 0:
                    P_b_cand = float(rb.P)
                    # Spike filter: in the retrograde region bubble P must not
                    # jump up; a >5 % rise above the last accepted value is a
                    # numerical artefact.
                    if bubble_pts and P_b_cand > bubble_pts[-1][1] * 1.05:
                        fine_bub_state = prev_bub_state
                    else:
                        P_b = P_b_cand
                        fine_bub_state = rb
                else:
                    fine_bub_state = None
            except Exception:
                fine_bub_state = None

            # ── dew ────────────────────────────────────────────────────────
            P_d: float | None = None
            try:
                rd = flash_obj.flash(zs=fracs, T=T_f, VF=1, hot_start=fine_dew_state)
                if rd.P is not None and rd.P > 0:
                    P_d = float(rd.P)
                    fine_dew_state = rd
                else:
                    fine_dew_state = None
            except Exception:
                fine_dew_state = None

            if P_b is not None and P_d is not None:
                frac_diff = abs(P_b - P_d) / max(P_b, P_d)
                if frac_diff < best_frac:
                    best_frac = frac_diff
                    best_close_T = T_f
                    best_close_P = (P_b + P_d) / 2.0
                if frac_diff < 0.05:
                    P_close_f = (P_b + P_d) / 2.0
                    _fill_cliff(bubble_pts, T_f, P_close_f)
                    bubble_pts.append((T_f, P_close_f))
                    dew_pts.append((T_f, P_close_f))
                    envelope_closed = True
                    break
                bubble_pts.append((T_f, P_b))
                dew_pts.append((T_f, P_d))
            elif P_b is not None:
                bubble_pts.append((T_f, P_b))
            elif P_d is not None:
                dew_pts.append((T_f, P_d))
            else:
                break  # both failed

            T_fine += fine_step

        # Force-close at the temperature of minimum divergence.
        if not envelope_closed and best_close_T is not None:
            bubble_pts = [(t, p) for t, p in bubble_pts if t < best_close_T]
            dew_pts    = [(t, p) for t, p in dew_pts    if t < best_close_T]

            # If the last bubble pressure is far above the closing pressure,
            # fill the gap with square-root interpolated points.  The scaling
            # P(T) = P_c + (P_last - P_c) × [(T_c - T)/(T_c - T_last)]^0.5
            # matches the mean-field critical exponent of cubic EOS and gives a
            # physically shaped closing nose rather than a vertical cliff.
            _fill_cliff(bubble_pts, best_close_T, best_close_P)

            bubble_pts.append((best_close_T, best_close_P))
            dew_pts.append((best_close_T, best_close_P))
        elif not envelope_closed and last_both is not None:
            T_fc2, P_fb, P_fd = last_both
            P_avg = (P_fb + P_fd) / 2.0
            bubble_pts = [(t, p) for t, p in bubble_pts if t <= T_fc2]
            dew_pts    = [(t, p) for t, p in dew_pts    if t <= T_fc2]
            if not bubble_pts or bubble_pts[-1][0] < T_fc2:
                bubble_pts.append((T_fc2, P_avg))
            if not dew_pts or dew_pts[-1][0] < T_fc2:
                dew_pts.append((T_fc2, P_avg))

    return bubble_pts, dew_pts, Tc, Pc
