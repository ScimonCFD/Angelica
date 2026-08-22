"""Phase envelope (bubble/dew point curves) computation for multi-component mixtures."""

from __future__ import annotations


def compute_phase_envelope(
    component_names: tuple[str, ...] | list[str],
    zs: tuple[float, ...] | list[float],
    n_T: int = 15,
) -> tuple[list[tuple[float, float]], list[tuple[float, float]], float, float]:
    """Compute bubble and dew point curves for a multi-component mixture.

    Args:
        component_names: Species identifiers accepted by thermo.Mixture.
        zs: Mole fractions (must sum to 1).
        n_T: Number of temperature points along each curve.

    Returns:
        (bubble_pts, dew_pts, Tc_K, Pc_Pa) where bubble_pts and dew_pts are
        lists of (T_K, P_Pa) ordered from low to high temperature, and
        Tc_K / Pc_Pa are the mole-fraction-weighted pseudo-critical point.
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
    T_hi = Tc * 0.97
    if T_lo >= T_hi:
        T_lo = T_hi * 0.50

    P_lo_search = 500.0
    P_hi_search = max(Pc * 3.0, 5.0e7)

    T_vals = np.linspace(T_lo, T_hi, n_T)

    def _vf(T: float, P: float) -> float | None:
        try:
            v = Mixture(names, zs=fracs, T=T, P=P).VF
            return float(v) if v is not None else None
        except Exception:
            return None

    def _bubble(T: float) -> float | None:
        v_hi = _vf(T, P_hi_search)
        if v_hi is None or v_hi > 0.01:
            return None
        v_lo = _vf(T, P_lo_search)
        if v_lo is None or v_lo < 0.01:
            return None
        lo, hi = P_lo_search, P_hi_search
        for _ in range(22):
            mid = (lo + hi) / 2.0
            v = _vf(T, mid)
            if v is None:
                hi = mid
            elif v > 0.01:
                lo = mid
            else:
                hi = mid
            if hi - lo < 200.0:
                break
        return (lo + hi) / 2.0

    def _dew(T: float) -> float | None:
        v_lo = _vf(T, P_lo_search)
        if v_lo is None or v_lo < 0.99:
            return None
        v_hi = _vf(T, P_hi_search)
        if v_hi is None or v_hi > 0.99:
            return None
        lo, hi = P_lo_search, P_hi_search
        for _ in range(22):
            mid = (lo + hi) / 2.0
            v = _vf(T, mid)
            if v is None:
                lo = mid
            elif v < 0.99:
                hi = mid
            else:
                lo = mid
            if hi - lo < 200.0:
                break
        return (lo + hi) / 2.0

    bubble_pts: list[tuple[float, float]] = []
    dew_pts: list[tuple[float, float]] = []

    for T in T_vals:
        P_b = _bubble(float(T))
        if P_b is not None:
            bubble_pts.append((float(T), P_b))
        P_d = _dew(float(T))
        if P_d is not None:
            dew_pts.append((float(T), P_d))

    return bubble_pts, dew_pts, Tc, Pc
