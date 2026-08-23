"""Phase envelope (bubble/dew point curves) computation for multi-component mixtures.

Uses thermo's built-in TVF flash (VF=0 bubble point, VF=1 dew point) which
applies the same Michelsen-style Newton iteration used by commercial simulators.
A closing bisection then finds the mixture critical point where the two curves meet.
"""

from __future__ import annotations


def compute_phase_envelope(
    component_names: tuple[str, ...] | list[str],
    zs: tuple[float, ...] | list[float],
    n_T: int = 20,
) -> tuple[list[tuple[float, float]], list[tuple[float, float]], float, float]:
    """Compute bubble and dew point curves for a multi-component mixture.

    Args:
        component_names: Species identifiers accepted by thermo.Mixture.
        zs: Mole fractions (must sum to 1).
        n_T: Number of temperature steps along each curve.

    Returns:
        (bubble_pts, dew_pts, Tc_K, Pc_Pa) where bubble_pts and dew_pts are
        lists of (T_K, P_Pa) ordered low→high temperature.  Both curves share
        the same last point (the mixture critical point) so the envelope is
        closed.  Tc_K / Pc_Pa are the mole-fraction-weighted pseudo-critical
        coordinates used only for the marker on the plot.
    """
    import numpy as np
    from thermo import Mixture

    names = list(component_names)
    fracs = list(zs)

    m0 = Mixture(names, zs=fracs, T=300.0, P=101325.0)
    Tcs: list[float] = list(m0.Tcs)
    Pcs: list[float] = list(m0.Pcs)
    Tbs: list[float] = list(m0.Tbs)

    Tc = sum(z * tc for z, tc in zip(fracs, Tcs))
    Pc = sum(z * pc for z, pc in zip(fracs, Pcs))

    T_lo = max(min(Tbs) * 0.65, 80.0)
    # Extend slightly past pseudo-Tc so the scan can reach convergence.
    T_hi_scan = Tc * 1.03
    if T_lo >= T_hi_scan:
        T_lo = T_hi_scan * 0.50

    P_lo_fallback = 500.0
    P_hi_fallback = max(Pc * 3.0, 5.0e7)

    # ── Flash helpers ──────────────────────────────────────────────────────────
    # Primary: ask thermo for the VF=0/VF=1 flash (bubble/dew point pressure).
    # Fallback: binary search on VF (works on older thermo builds).

    def _bubble(T: float) -> float | None:
        T = float(T)
        try:
            m = Mixture(names, zs=fracs, T=T, VF=0)
            if m.P is not None and m.P > 0:
                return float(m.P)
        except Exception:
            pass
        # Fallback: binary search
        return _binary_bubble(T)

    def _dew(T: float) -> float | None:
        T = float(T)
        try:
            m = Mixture(names, zs=fracs, T=T, VF=1)
            if m.P is not None and m.P > 0:
                return float(m.P)
        except Exception:
            pass
        # Fallback: binary search
        return _binary_dew(T)

    def _vf(T: float, P: float) -> float | None:
        try:
            v = Mixture(names, zs=fracs, T=T, P=P).VF
            return float(v) if v is not None else None
        except Exception:
            return None

    def _binary_bubble(T: float) -> float | None:
        v_hi = _vf(T, P_hi_fallback)
        if v_hi is None or v_hi > 0.01:
            return None
        v_lo = _vf(T, P_lo_fallback)
        if v_lo is None or v_lo < 0.01:
            return None
        lo, hi = P_lo_fallback, P_hi_fallback
        for _ in range(24):
            mid = (lo + hi) / 2.0
            v = _vf(T, mid)
            if v is None:
                hi = mid
            elif v > 0.01:
                lo = mid
            else:
                hi = mid
            if hi - lo < 100.0:
                break
        return (lo + hi) / 2.0

    def _binary_dew(T: float) -> float | None:
        v_lo = _vf(T, P_lo_fallback)
        if v_lo is None or v_lo < 0.99:
            return None
        v_hi = _vf(T, P_hi_fallback)
        if v_hi is None or v_hi > 0.99:
            return None
        lo, hi = P_lo_fallback, P_hi_fallback
        for _ in range(24):
            mid = (lo + hi) / 2.0
            v = _vf(T, mid)
            if v is None:
                lo = mid
            elif v < 0.99:
                hi = mid
            else:
                lo = mid
            if hi - lo < 100.0:
                break
        return (lo + hi) / 2.0

    # ── Temperature scan ───────────────────────────────────────────────────────
    T_vals = np.linspace(T_lo, T_hi_scan, n_T + 8)

    bubble_pts: list[tuple[float, float]] = []
    dew_pts: list[tuple[float, float]] = []
    envelope_closed = False

    for T in T_vals:
        T = float(T)
        P_b = _bubble(T)
        P_d = _dew(T)

        if P_b is not None and P_d is not None:
            frac_diff = abs(P_b - P_d) / max(P_b, P_d)
            if frac_diff < 0.03:
                # Curves have converged — critical point reached.
                P_crit_approx = (P_b + P_d) / 2.0
                bubble_pts.append((T, P_crit_approx))
                dew_pts.append((T, P_crit_approx))
                envelope_closed = True
                break
            bubble_pts.append((T, P_b))
            dew_pts.append((T, P_d))
        elif P_b is not None:
            bubble_pts.append((T, P_b))
        elif P_d is not None:
            dew_pts.append((T, P_d))
        elif bubble_pts or dew_pts:
            # Both searches failed — past the two-phase region.
            break

    # ── Close the envelope (bisection toward the critical point) ──────────────
    # If the scan ended without the two curves meeting, bisect in the gap
    # above the last scan temperature to find where P_bubble ≈ P_dew.
    if not envelope_closed and bubble_pts and dew_pts:
        T_cl_lo = max(bubble_pts[-1][0], dew_pts[-1][0])
        T_cl_hi = Tc * 1.05

        for _ in range(16):
            if T_cl_hi - T_cl_lo < 0.02:
                break
            T_mid = (T_cl_lo + T_cl_hi) / 2.0
            P_b = _bubble(T_mid)
            P_d = _dew(T_mid)

            if P_b is None or P_d is None:
                T_cl_hi = T_mid
                continue

            frac_diff = abs(P_b - P_d) / max(P_b, P_d)
            if frac_diff < 0.03:
                P_crit_approx = (P_b + P_d) / 2.0
                bubble_pts = [(t, p) for t, p in bubble_pts if t < T_mid]
                dew_pts = [(t, p) for t, p in dew_pts if t < T_mid]
                bubble_pts.append((T_mid, P_crit_approx))
                dew_pts.append((T_mid, P_crit_approx))
                break
            T_cl_lo = T_mid  # curves still separated — move up

    return bubble_pts, dew_pts, Tc, Pc
