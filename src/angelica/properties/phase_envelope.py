"""Phase envelope (bubble/dew point curves) via Michelsen (1980) arc-length.

Algorithm: Wilson K-value bootstrap finds the first two-phase point without
scanning.  Arc-length continuation (predictor = SVD null vector; corrector =
Newton-Raphson) then traces the full envelope in both directions from that
starting point.  Falls back to a temperature scan when Wilson fails.
"""

from __future__ import annotations

from typing import Any

import numpy as np


def _build_kij_matrix(constants) -> list[list[float]]:
    """Build kij matrix from EPPR78 database (Jaubert & Mutelet 2004).

    Returns an N×N list of lists; off-diagonal values come from thermo's
    bundled eppr78_common.json keyed by CAS-number pairs.  Missing pairs
    default to 0.0.
    """
    import json
    import os

    try:
        import thermo as _thermo_pkg
        db_path = os.path.join(
            os.path.dirname(_thermo_pkg.__file__),
            "Interaction Parameters",
            "eppr78_common.json",
        )
        with open(db_path, encoding="utf-8") as fh:
            raw: dict = json.load(fh)["data"]
    except Exception:
        raw = {}

    CASs = list(constants.CASs)
    N = len(CASs)
    kijs: list[list[float]] = [[0.0] * N for _ in range(N)]
    for i in range(N):
        for j in range(i + 1, N):
            pair = raw.get(f"{CASs[i]} {CASs[j]}") or raw.get(f"{CASs[j]} {CASs[i]}")
            kij = float(pair["kij"]) if pair else 0.0
            kijs[i][j] = kij
            kijs[j][i] = kij
    return kijs


def _wilson_start(
    flash_obj: Any,
    constants: Any,
    z: np.ndarray,
    fracs: list[float],
    T_lo: float,
    T_hi: float,
    VF: int,
) -> tuple[np.ndarray, float, float] | None:
    """Find a valid two-phase starting point using the Wilson K-value correlation.

    Wilson (1968) gives an analytical estimate of bubble/dew pressure:
      K_i = (Pci/P) * exp(5.373*(1+ωi)*(1 - Tci/T))
      Bubble: P_bub = Σ zi*Pci*exp(...)
      Dew:    1/P_dew = Σ zi/(Pci*exp(...))

    A single TP flash at the Wilson (T, P) estimate converges immediately
    because the composition is already near the phase boundary.  This
    replaces the brute-force temperature scan (~200-400 flash calls) with
    1-5 flash calls.

    Returns (K, T, P) for a valid two-phase point, or None if all attempts
    fail (caller falls back to the legacy scan).
    """
    Tcs_ = np.array(constants.Tcs)
    Pcs_ = np.array(constants.Pcs)
    omegas_ = np.array(constants.omegas)
    z_ = np.array(fracs)
    N = len(fracs)
    Tc_mix = float(np.dot(z_, Tcs_))

    # Candidate fractions of Tc_mix: bubble curve is well-conditioned
    # around 0.6-0.7*Tc; dew curve first appears at lower T (left leg).
    if VF == 0:
        candidates = (0.65, 0.55, 0.75, 0.50, 0.80)
    else:
        candidates = (0.50, 0.45, 0.55, 0.40, 0.60)

    for frac_T in candidates:
        T_try = frac_T * Tc_mix
        if not (T_lo + 1.0 < T_try < T_hi - 1.0):
            continue

        exp_term = np.exp(5.373 * (1.0 + omegas_) * (1.0 - Tcs_ / T_try))

        if VF == 0:
            P_w = float(np.dot(z_, Pcs_ * exp_term))
        else:
            denom = Pcs_ * exp_term
            with np.errstate(divide="ignore", invalid="ignore"):
                inv_sum = float(np.dot(z_, np.where(denom > 0, 1.0 / denom, 0.0)))
            P_w = 1.0 / max(inv_sum, 1e-30)

        if not (1e3 < P_w < 1e8):
            continue

        try:
            res = flash_obj.flash(zs=fracs, T=float(T_try), P=float(P_w))
            if res is None or res.phase_count != 2 or not res.P or res.P <= 1e3:
                continue
            K = np.array([res.gas.zs[i] / res.liquid0.zs[i] for i in range(N)])
            if not np.all(np.isfinite(K)) or not np.all(K > 0):
                continue
            if float(np.max(np.abs(np.log(K)))) > 0.3:
                return K, float(T_try), float(res.P)
        except Exception:
            pass

    return None


def _eval_FJ(
    gas_phase: Any,
    liq_phase: Any,
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
                if np.all(np.isfinite(F_t)) and np.max(np.abs(F_t)) < err * (1 - 1e-6):
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
    ascend: bool = True,
    T_min: float = 50.0,
    P_max: float = 1e8,
) -> "list[tuple[float, float]]":
    """Arc-length continuation from (T_start, P_start, K_start).

    ascend=True  traces toward higher T (toward the critical point).
    ascend=False traces toward lower T (away from the critical point,
                 toward the left leg of the envelope).  Stops when
                 T < T_min or P > P_max.
    """
    N = len(z)
    S = np.concatenate(
        [np.log(np.clip(K_start, 1e-15, 1e15)), [np.log(T_start), np.log(P_start)]]
    )

    # Initial tangent seed: +T direction for ascent, −T for descent.
    tangent = np.zeros(N + 2)
    tangent[N] = 1.0 if ascend else -1.0
    try:
        _, J0 = _eval_FJ(gas_phase, liq_phase, T_start, P_start, K_start, z, VF)
        tangent = _null(J0, tangent)
        if ascend:
            if tangent[N] < 0:
                tangent = -tangent
        else:
            if tangent[N] > 0:
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

        # Ascent: stop near the critical point (K → 1).
        if ascend and lnK_max_new < lnK_stop:
            break

        # Descent: stop when we exit the practical T/P range.
        if not ascend and (T_new < T_min or P_new > P_max):
            break

    return pts


def _bubble_trace(
    gas_phase: Any,
    liq_phase: Any,
    flash_obj: Any,
    constants: Any,
    z: "np.ndarray",
    fracs: list[float],
    T_lo: float,
    T_hi: float,
    lnK_arc_start: float = 0.5,
) -> "list[tuple[float, float]]":
    """Bubble-point curve: Wilson bootstrap → bidirectional arc-length."""
    N = len(z)

    # ── Wilson bootstrap ──────────────────────────────────────────────────
    ws = _wilson_start(flash_obj, constants, z, fracs, T_lo, T_hi, VF=0)
    if ws is not None:
        K_start, T_start, P_start = ws
        if float(np.max(np.abs(np.log(K_start)))) >= lnK_arc_start:
            pts_up   = _trace_arc(gas_phase, liq_phase, z, 0, T_start, P_start, K_start,
                                   ascend=True)
            pts_down = _trace_arc(gas_phase, liq_phase, z, 0, T_start, P_start, K_start,
                                   ascend=False, T_min=T_lo)
            # pts_down is [T_start, T_lower_1, ...]: reverse to get low→high T order.
            return list(reversed(pts_down)) + pts_up[1:]

    # ── Fallback: brute-force temperature scan ────────────────────────────
    prev = None
    pts: list[tuple[float, float]] = []
    arc_K: "np.ndarray | None" = None
    arc_P = arc_T = 0.0

    for T in np.arange(T_lo, T_hi, 0.5):
        try:
            hs: dict[str, Any] = {"hot_start": prev} if prev is not None else {}
            res = flash_obj.flash(zs=fracs, T=float(T), VF=0, **hs)
            if res.phase_count != 2 or not res.P or res.P <= 1e3:
                continue
            K = np.array([res.gas.zs[i] / res.liquid0.zs[i] for i in range(N)])
            if not np.all(np.isfinite(K)) or not np.all(K > 0):
                continue
            lnK_max = float(np.max(np.abs(np.log(K))))

            prev_lkm = float(np.max(np.abs(np.log(arc_K)))) if arc_K is not None else 1.0
            if arc_K is not None and lnK_max < 0.01 and prev_lkm > 0.1:
                break

            if pts and abs(res.P - pts[-1][1]) > pts[-1][1] * 0.5:
                break

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
    gas_phase: Any,
    liq_phase: Any,
    flash_obj: Any,
    constants: Any,
    z: "np.ndarray",
    fracs: list[float],
    T_lo: float,
    T_hi: float,
) -> "list[tuple[float, float]]":
    """Dew-point curve: Wilson bootstrap → bidirectional arc-length."""
    N = len(z)

    # ── Wilson bootstrap ──────────────────────────────────────────────────
    ws = _wilson_start(flash_obj, constants, z, fracs, T_lo, T_hi, VF=1)
    if ws is not None:
        K_start, T_start, P_start = ws
        if float(np.max(np.abs(np.log(K_start)))) > 0.3:
            pts_up   = _trace_arc(gas_phase, liq_phase, z, 1, T_start, P_start, K_start,
                                   ascend=True)
            pts_down = _trace_arc(gas_phase, liq_phase, z, 1, T_start, P_start, K_start,
                                   ascend=False, T_min=T_lo)
            return list(reversed(pts_down)) + pts_up[1:]

    # ── Fallback: scan ────────────────────────────────────────────────────
    prev = None
    pts: list[tuple[float, float]] = []
    arc_K: "np.ndarray | None" = None
    arc_P = arc_T = 0.0

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
                break

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
    n_T: int = 20,  # noqa: ARG001  accepted for API compatibility, arc-length is used instead
    eos_name: str = "PR",
) -> "tuple[list[tuple[float, float]], list[tuple[float, float]], float, float]":
    """Compute bubble and dew point curves for a multi-component mixture.

    Returns:
        (bubble_pts, dew_pts, Tc_K, Pc_Pa) ordered low→high temperature.
        Both curves share the same last point (the critical point) so the
        envelope is closed.  Tc_K / Pc_Pa are mole-fraction-weighted
        pseudo-critical coordinates used only for the plot marker.
        n_T is accepted for API compatibility but ignored (arc-length is used).
        eos_name: ``"PR"`` (Peng-Robinson, default) or ``"SRK"``.
    """
    import warnings

    from thermo import ChemicalConstantsPackage
    from thermo.eos_mix import PRMIX, SRKMIX
    from thermo.flash import FlashVL
    from thermo.phases import CEOSGas, CEOSLiquid

    eos_cls = SRKMIX if eos_name.upper() == "SRK" else PRMIX

    names = list(component_names)
    fracs = list(zs)
    z = np.asarray(fracs, float)

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        constants, props = ChemicalConstantsPackage.from_IDs(names)
        kijs = _build_kij_matrix(constants)
        eos_kw = {"Tcs": constants.Tcs, "Pcs": constants.Pcs, "omegas": constants.omegas, "kijs": kijs}
        gas_phase = CEOSGas(eos_cls, eos_kw, HeatCapacityGases=props.HeatCapacityGases)
        liq_phase = CEOSLiquid(eos_cls, eos_kw, HeatCapacityGases=props.HeatCapacityGases)
        flash_obj: Any = FlashVL(constants, props, liquid=liq_phase, gas=gas_phase)

    Tcs: list[float] = list(constants.Tcs)
    Pcs: list[float] = list(constants.Pcs)
    Tbs: list[float] = list(constants.Tbs)

    Tc = sum(zi * tc for zi, tc in zip(fracs, Tcs))
    Pc = sum(zi * pc for zi, pc in zip(fracs, Pcs))

    T_lo = max(min(Tbs) * 0.55, 75.0)
    T_hi = Tc * 1.5

    bubble_pts = _bubble_trace(gas_phase, liq_phase, flash_obj, constants, z, fracs, T_lo, T_hi)
    dew_pts = _dew_trace(gas_phase, liq_phase, flash_obj, constants, z, fracs, T_lo, T_hi)

    if bubble_pts and dew_pts:
        T_bc, P_bc = bubble_pts[-1]
        T_dc, P_dc = dew_pts[-1]
        T_crit = (T_bc + T_dc) / 2.0
        P_crit = (P_bc + P_dc) / 2.0
        Tc, Pc = T_crit, P_crit
        if abs(T_bc - T_dc) < 2.0 and abs(P_bc - P_dc) / max(P_bc, P_dc) < 0.02:
            bubble_pts[-1] = (T_crit, P_crit)
            dew_pts[-1] = (T_crit, P_crit)

    return bubble_pts, dew_pts, Tc, Pc


def compute_quality_line(
    component_names: "tuple[str, ...] | list[str]",
    zs: "tuple[float, ...] | list[float]",
    vf: float,
    eos_name: str = "PR",
) -> "list[tuple[float, float]]":
    """Compute a constant-vapor-fraction (quality) line inside the phase envelope.

    Returns [(T_K, P_Pa)] ordered from low to high temperature, covering the
    two-phase region at the given vapor fraction (0 < vf < 1).
    eos_name: ``"PR"`` (Peng-Robinson, default) or ``"SRK"``.
    """
    import warnings

    from thermo import ChemicalConstantsPackage
    from thermo.eos_mix import PRMIX, SRKMIX
    from thermo.flash import FlashVL
    from thermo.phases import CEOSGas, CEOSLiquid

    eos_cls = SRKMIX if eos_name.upper() == "SRK" else PRMIX

    names = list(component_names)
    fracs = list(zs)

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        constants, props = ChemicalConstantsPackage.from_IDs(names)
        kijs = _build_kij_matrix(constants)
        eos_kw = {"Tcs": constants.Tcs, "Pcs": constants.Pcs, "omegas": constants.omegas, "kijs": kijs}
        gas_phase = CEOSGas(eos_cls, eos_kw, HeatCapacityGases=props.HeatCapacityGases)
        liq_phase = CEOSLiquid(eos_cls, eos_kw, HeatCapacityGases=props.HeatCapacityGases)
        flash_obj: Any = FlashVL(constants, props, liquid=liq_phase, gas=gas_phase)

    Tcs_list: list[float] = list(constants.Tcs)
    Tbs: list[float] = list(constants.Tbs)
    Tc_mix = sum(zi * tc for zi, tc in zip(fracs, Tcs_list))
    T_lo = max(min(Tbs) * 0.55, 75.0)
    T_hi = Tc_mix * 1.5

    # Use Wilson to estimate a good starting temperature: the quality line
    # starts near the bubble curve (for low vf) or dew curve (for high vf).
    # A fraction of Tc between 0.4 and 0.6 is reliably inside the envelope.
    T_scan_start = max(T_lo, min(0.5 * Tc_mix, T_hi - 5.0))

    pts: list[tuple[float, float]] = []
    prev = None

    # Scan with 2 °C steps (vs. original 1 °C) starting from Wilson estimate.
    for T in np.arange(T_scan_start, T_hi, 2.0):
        try:
            hs: dict[str, Any] = {"hot_start": prev} if prev is not None else {}
            res = flash_obj.flash(zs=fracs, T=float(T), VF=vf, **hs)
            if res.phase_count == 2 and res.P and res.P > 1e3:
                pts.append((float(T), float(res.P)))
                prev = res
            elif pts:
                break
        except Exception:
            pass

    # If nothing found from T_scan_start, fall back to scanning from T_lo.
    if not pts:
        prev = None
        for T in np.arange(T_lo, T_scan_start, 2.0):
            try:
                hs_: dict[str, Any] = {"hot_start": prev} if prev is not None else {}
                res = flash_obj.flash(zs=fracs, T=float(T), VF=vf, **hs_)
                if res.phase_count == 2 and res.P and res.P > 1e3:
                    pts.append((float(T), float(res.P)))
                    prev = res
                elif pts:
                    break
            except Exception:
                pass

    return pts
