"""Phase envelope (bubble/dew point curves) via Michelsen (1980) arc-length.

Algorithm: hot-start flash scan up to the near-critical region, then arc-length
continuation (predictor = SVD null vector; corrector = Newton-Raphson) to trace
through the critical point.  Works for both bubble (VF=0) and dew (VF=1).
"""

from __future__ import annotations

import numpy as np


def _eval_FJ(
    gas_phase: object,
    liq_phase: object,
    T: float,
    P: float,
    K: "np.ndarray",
    z: "np.ndarray",
    VF: int,
) -> "tuple[np.ndarray, np.ndarray]":
    """Michelsen equilibrium equations F and their Jacobian J.

    State vector S = [lnK_1..lnK_N, lnT, lnP] (length N+2).
    Equations  F = [lnK_i - (lnφL_i - lnφV_i) for all i, ln(S_)] (length N+1).
    Jacobian   J has shape (N+1, N+2).

    For bubble (VF=0): x=z, y=Kz/S_, S_=Σ K_i z_i
    For dew   (VF=1): y=z, x=z/(K S_), S_=Σ z_i/K_i
    """
    N = len(z)
    if VF == 0:
        S_ = float(np.dot(K, z))
        x, y = z, K * z / S_
    else:
        S_ = float(np.dot(z / K, np.ones(N)))
        y, x = z, z / (K * S_)

    liq_obj = liq_phase.to(T=T, P=P, zs=list(x))
    gas_obj = gas_phase.to(T=T, P=P, zs=list(y))

    lnφL = np.array(liq_obj.lnphis())
    lnφV = np.array(gas_obj.lnphis())
    F = np.empty(N + 1)
    F[:N] = np.log(K) - (lnφL - lnφV)
    F[N] = np.log(abs(S_))

    dlnφL_dT = np.array(liq_obj.dlnphis_dT())
    dlnφV_dT = np.array(gas_obj.dlnphis_dT())
    dlnφL_dP = np.array(liq_obj.dlnphis_dP())
    dlnφV_dP = np.array(gas_obj.dlnphis_dP())

    if VF == 0:
        D, q, sgn = np.array(gas_obj.dlnphis_dzs()), y, +1.0
    else:
        D, q, sgn = np.array(liq_obj.dlnphis_dzs()), x, -1.0

    A = D * q[np.newaxis, :]
    C = A.sum(axis=1)

    J = np.zeros((N + 1, N + 2))
    np.fill_diagonal(J[:N, :N], 1.0)
    J[:N, :N] += A
    J[:N, :N] -= np.outer(C, q)
    J[:N, N]   = -T * (dlnφL_dT - dlnφV_dT)
    J[:N, N+1] = -P * (dlnφL_dP - dlnφV_dP)
    J[N,  :N]  = sgn * q
    return F, J


def _null(J: "np.ndarray", prev_t: "np.ndarray") -> "np.ndarray":
    """Return the right null vector of J consistent with prev_t."""
    _, _, Vt = np.linalg.svd(J, full_matrices=True)
    t = Vt[-1, :].copy()
    if np.dot(t, prev_t) < 0:
        t = -t
    norm = np.linalg.norm(t)
    return t / (norm if norm > 1e-300 else 1.0)


def _correct(
    gas_phase: object,
    liq_phase: object,
    T0: float,
    P0: float,
    K0: "np.ndarray",
    z: "np.ndarray",
    VF: int,
    spec_idx: int,
    spec_val: float,
    tol: float = 1e-9,
    max_iter: int = 50,
) -> "np.ndarray | None":
    """Newton-Raphson corrector with line search.  Returns S or None on failure."""
    N = len(z)
    S = np.concatenate(
        [
            np.log(np.clip(K0, 1e-15, 1e15)),
            [np.log(max(T0, 50.0)), np.log(max(P0, 1e3))],
        ]
    )
    S[spec_idx] = spec_val

    err = float("inf")
    for _ in range(max_iter):
        K_ = np.exp(np.clip(S[:N], -40, 40))
        T_ = max(float(np.exp(np.clip(S[N], 4.5, 8.0))), 90.0)
        P_ = max(float(np.exp(np.clip(S[N + 1], 7.0, 19.0))), 1e3)
        try:
            F, J = _eval_FJ(gas_phase, liq_phase, T_, P_, K_, z, VF)
        except Exception:
            return None
        if not np.all(np.isfinite(F)):
            return None
        err = float(np.max(np.abs(F)))
        if err < tol:
            return S.copy()

        F_aug = np.append(F, S[spec_idx] - spec_val)
        J_aug = np.zeros((N + 2, N + 2))
        J_aug[:N + 1, :] = J
        J_aug[N + 1, spec_idx] = 1.0
        try:
            dS = np.linalg.solve(J_aug, -F_aug)
        except np.linalg.LinAlgError:
            return None
        if not np.all(np.isfinite(dS)):
            return None

        alpha = 1.0
        for _ in range(14):
            S_try = S + alpha * dS
            S_try[spec_idx] = spec_val
            K_t = np.exp(np.clip(S_try[:N], -40, 40))
            T_t = max(float(np.exp(np.clip(S_try[N], 4.5, 8.0))), 90.0)
            P_t = max(float(np.exp(np.clip(S_try[N + 1], 7.0, 19.0))), 1e3)
            try:
                F_t, _ = _eval_FJ(gas_phase, liq_phase, T_t, P_t, K_t, z, VF)
                if np.all(np.isfinite(F_t)) and np.max(np.abs(F_t)) < err * (1 - 1e-4):
                    S = S_try
                    break
            except Exception:
                pass
            alpha *= 0.5
        else:
            return None

    return S if err < 1e-5 else None


def _trace_arc(
    gas_phase: object,
    liq_phase: object,
    z: "np.ndarray",
    VF: int,
    T_start: float,
    P_start: float,
    K_start: "np.ndarray",
    ds0: float = 0.05,
    ds_min: float = 0.001,
    ds_max: float = 0.5,
    max_steps: int = 1000,
    lnK_stop: float = 0.005,
) -> "list[tuple[float, float]]":
    """Arc-length continuation from (T_start, P_start, K_start) toward the critical."""
    N = len(z)
    S = np.concatenate(
        [np.log(np.clip(K_start, 1e-15, 1e15)), [np.log(T_start), np.log(P_start)]]
    )

    tangent = np.zeros(N + 2)
    tangent[N] = 1.0
    try:
        _, J0 = _eval_FJ(gas_phase, liq_phase, T_start, P_start, K_start, z, VF)
        tangent = _null(J0, tangent)
        if tangent[N] < 0:
            tangent = -tangent
    except Exception:
        pass

    pts: list[tuple[float, float]] = [(T_start, P_start)]
    ds = ds0
    n_fail = 0

    for _ in range(max_steps):
        K_cur = np.exp(np.clip(S[:N], -40, 40))
        lnK_max = float(np.max(np.abs(np.log(K_cur))))

        ds = min(ds, lnK_max * 0.15, ds_max)
        ds = max(ds, ds_min)

        S_pred = S + ds * tangent
        spec_idx = int(np.argmax(np.abs(tangent)))
        spec_val = float(S_pred[spec_idx])

        K_pred = np.exp(np.clip(S_pred[:N], -40, 40))
        T_pred = max(float(np.exp(np.clip(S_pred[N], 4.5, 8.0))), 90.0)
        P_pred = max(float(np.exp(np.clip(S_pred[N + 1], 7.0, 19.0))), 1e3)

        S_new = _correct(
            gas_phase, liq_phase, T_pred, P_pred, K_pred, z, VF, spec_idx, spec_val
        )

        if S_new is None:
            ds = max(ds * 0.5, ds_min)
            n_fail += 1
            if n_fail > 60:
                break
            continue

        n_fail = 0
        K_new = np.exp(np.clip(S_new[:N], -40, 40))
        T_new = float(np.exp(np.clip(S_new[N], 4.5, 8.0)))
        P_new = float(np.exp(np.clip(S_new[N + 1], 7.0, 19.0)))
        lnK_max_new = float(np.max(np.abs(np.log(K_new))))

        try:
            _, Jn = _eval_FJ(gas_phase, liq_phase, T_new, P_new, K_new, z, VF)
            tangent = _null(Jn, tangent)
        except Exception:
            pass

        S = S_new
        pts.append((T_new, P_new))

        if lnK_max_new > 0.3:
            ds = min(ds * 1.05, ds_max)

        if lnK_max_new < lnK_stop:
            break

    return pts


def _bubble_trace(
    gas_phase: object,
    liq_phase: object,
    flash_obj: object,
    z: "np.ndarray",
    names: list[str],
    fracs: list[float],
    T_lo: float,
    T_hi: float,
    lnK_arc_start: float = 0.5,
) -> "list[tuple[float, float]]":
    """Bubble-point curve via flash scan then arc-length continuation."""
    N = len(z)
    prev = None
    pts: list[tuple[float, float]] = []
    arc_K: "np.ndarray | None" = None
    arc_P = arc_T = 0.0

    for T in np.arange(T_lo, T_hi, 0.5):
        try:
            hs = {"hot_start": prev} if prev is not None else {}
            res = flash_obj.flash(zs=fracs, T=float(T), VF=0, **hs)
            if res.phase_count != 2 or not res.P or res.P <= 1e3:
                continue
            K = np.array([res.gas.zs[i] / res.liquid0.zs[i] for i in range(N)])
            if not np.all(np.isfinite(K)) or not np.all(K > 0):
                continue
            lnK_max = float(np.max(np.abs(np.log(K))))

            prev_lkm = float(np.max(np.abs(np.log(arc_K)))) if arc_K is not None else 1.0
            if arc_K is not None and lnK_max < 0.01 and prev_lkm > 0.1:
                break  # snap to trivial — switch to arc-length

            if pts and abs(res.P - pts[-1][1]) > pts[-1][1] * 0.5:
                break  # pressure jump — switch to arc-length

            pts.append((float(T), float(res.P)))
            prev = res
            if lnK_max >= lnK_arc_start:
                arc_K = K.copy()
                arc_P = float(res.P)
                arc_T = float(T)
        except Exception:
            pass

    if arc_K is not None:
        pts = [(t, p) for t, p in pts if t <= arc_T]
        arc_pts = _trace_arc(gas_phase, liq_phase, z, 0, arc_T, arc_P, arc_K)
        pts.extend(arc_pts[1:])

    return pts


def _dew_trace(
    gas_phase: object,
    liq_phase: object,
    flash_obj: object,
    z: "np.ndarray",
    names: list[str],
    fracs: list[float],
    T_lo: float,
    T_hi: float,
) -> "list[tuple[float, float]]":
    """Dew-point curve via flash scan then arc-length continuation."""
    N = len(z)
    prev = None
    pts: list[tuple[float, float]] = []
    arc_K: "np.ndarray | None" = None
    arc_P = arc_T = 0.0

    # find a well-conditioned starting point
    for T_s in np.arange(T_lo, T_lo + 100.0, 5.0):
        try:
            res = flash_obj.flash(zs=fracs, T=float(T_s), VF=1)
            if res.phase_count == 2 and res.P and res.P > 1e3:
                K = np.array([res.gas.zs[i] / res.liquid0.zs[i] for i in range(N)])
                if (
                    np.all(np.isfinite(K))
                    and np.all(K > 0)
                    and float(np.max(np.abs(np.log(K)))) > 0.5
                ):
                    prev = res
                    arc_K = K.copy()
                    arc_P = float(res.P)
                    arc_T = float(T_s)
                    pts.append((float(T_s), float(res.P)))
                    break
        except Exception:
            pass

    if arc_K is None:
        return []

    for T in np.arange(arc_T + 0.5, T_hi, 0.5):
        try:
            res = flash_obj.flash(zs=fracs, T=float(T), VF=1, hot_start=prev)
            if res.phase_count != 2 or not res.P or res.P <= 1e3:
                continue
            K = np.array([res.gas.zs[i] / res.liquid0.zs[i] for i in range(N)])
            if not np.all(np.isfinite(K)) or not np.all(K > 0):
                continue
            lnK_max = float(np.max(np.abs(np.log(K))))

            prev_lkm = float(np.max(np.abs(np.log(arc_K)))) if arc_K is not None else 1.0
            if arc_K is not None and lnK_max < 0.01 and prev_lkm > 0.1:
                break  # snap — switch to arc-length

            pts.append((float(T), float(res.P)))
            arc_K = K.copy()
            arc_P = float(res.P)
            arc_T = float(T)
            prev = res
        except Exception:
            pass

    arc_pts = _trace_arc(gas_phase, liq_phase, z, 1, arc_T, arc_P, arc_K)
    pts.extend(arc_pts[1:])

    return pts


def compute_phase_envelope(
    component_names: "tuple[str, ...] | list[str]",
    zs: "tuple[float, ...] | list[float]",
    n_T: int = 20,
) -> "tuple[list[tuple[float, float]], list[tuple[float, float]], float, float]":
    """Compute bubble and dew point curves for a multi-component mixture.

    Returns:
        (bubble_pts, dew_pts, Tc_K, Pc_Pa) ordered low→high temperature.
        Both curves share the same last point (the critical point) so the
        envelope is closed.  Tc_K / Pc_Pa are mole-fraction-weighted
        pseudo-critical coordinates used only for the plot marker.
        n_T is accepted for API compatibility but ignored (arc-length is used).
    """
    import warnings

    warnings.filterwarnings("ignore")

    from thermo import ChemicalConstantsPackage
    from thermo.flash import FlashVL
    from thermo.phases import CEOSGas, CEOSLiquid
    from thermo.eos_mix import PRMIX

    names = list(component_names)
    fracs = list(zs)
    z = np.asarray(fracs, float)

    constants, props = ChemicalConstantsPackage.from_IDs(names)
    eos_kw = dict(Tcs=constants.Tcs, Pcs=constants.Pcs, omegas=constants.omegas)
    gas_phase = CEOSGas(PRMIX, eos_kw, HeatCapacityGases=props.HeatCapacityGases)
    liq_phase = CEOSLiquid(PRMIX, eos_kw, HeatCapacityGases=props.HeatCapacityGases)
    flash_obj = FlashVL(constants, props, liquid=liq_phase, gas=gas_phase)

    Tcs: list[float] = list(constants.Tcs)
    Pcs: list[float] = list(constants.Pcs)
    Tbs: list[float] = list(constants.Tbs)

    Tc = sum(zi * tc for zi, tc in zip(fracs, Tcs))
    Pc = sum(zi * pc for zi, pc in zip(fracs, Pcs))

    T_lo = max(min(Tbs) * 0.55, 75.0)
    T_hi = Tc * 1.5

    bubble_pts = _bubble_trace(
        gas_phase, liq_phase, flash_obj, z, names, fracs, T_lo, T_hi
    )
    dew_pts = _dew_trace(
        gas_phase, liq_phase, flash_obj, z, names, fracs, T_lo, T_hi
    )

    # Close the envelope at the critical: use the last point of each arc that
    # is closest to the other.  If both arcs reach lnK<0.005 they've converged;
    # average the final T/P pair as the shared critical point.
    if bubble_pts and dew_pts:
        T_bc, P_bc = bubble_pts[-1]
        T_dc, P_dc = dew_pts[-1]
        T_crit = (T_bc + T_dc) / 2.0
        P_crit = (P_bc + P_dc) / 2.0
        # Both arcs should be within 1 K / 1 bar — if so, merge; otherwise keep as-is
        if abs(T_bc - T_dc) < 2.0 and abs(P_bc - P_dc) / max(P_bc, P_dc) < 0.02:
            bubble_pts[-1] = (T_crit, P_crit)
            dew_pts[-1] = (T_crit, P_crit)

    return bubble_pts, dew_pts, Tc, Pc
