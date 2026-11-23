"""
Gas-hydrate formation conditions — van der Waals–Platteeuw model.

Equilibrium of water between liquid water (L) and hydrate (H), with the
empty hydrate lattice (β) as intermediate state:

    Δμ_w^{β−L}(T, P) = Δμ_w^{β−H}(T, P)

Hydrate side (van der Waals & Platteeuw 1959):

    Δμ_w^{β−H} / RT = −Σ_m ν_m · ln(1 − Σ_j θ_mj),
    θ_mj = C_mj·f_j / (1 + Σ_k C_mk·f_k),     C_mj = (A_mj / T)·exp(B_mj / T)

with guest fugacities f_j (atm) from the Peng-Robinson EOS of the hydrate
formers (dry basis; water is assumed present in excess).  ν: sI small 1/23,
large 3/23; sII small 2/17, large 1/17.

Water side (Holder et al. 1980 form), liquid-water basis for T ≥ T0:

    Δμ_w^{β−L} / RT = Δμ0 / (R·T0) − ∫_{T0}^{T} Δh_w / (R·T²) dT + Δv_w·P / (R·T) − ln x_w
    Δh_w(T) = Δh0_L + ∫_{T0}^{T} Δcp dT,   Δcp = −38.12 + 0.141·(T − T0)  J/(mol·K)

Parameters
----------
Langmuir constants: Munck, Skjold-Jørgensen & Rasmussen (1988), Chem. Eng.
Sci. 43, 2661, as tabulated (A in K/atm, B in K) by Antunes et al. (2018),
Braz. J. Chem. Eng. 35, 265, Table 2.

Empty-lattice reference properties (T0 = 273.15 K): sI Δμ0 = 1264 J/mol,
Δh0 = 1151 J/mol (ice), Δv0 = 3.0 cm³/mol (ice); sII 883 J/mol, 808 J/mol,
3.4 cm³/mol (Holder et al. 1980, the set Munck's constants were fitted with;
values as tabulated by Dongre et al. 2022, Sci. Rep. 12, 16688, Table 1).
Liquid basis: Δh0_L = Δh0 − 6009.5 J/mol (ice fusion), Δv0_L = Δv0 + 1.6 cm³/mol.
Δcp from Holder & Hand (1982), as given by Dongre et al. (2022), Eq. 7.

Validation status
-----------------
Quadruple points of CH4, C2H6, C3H8 and CO2 (Sloan & Koh values seen through
a secondary source): within 7.3 % in pressure at Q1 and 0.8 K at Q2.  Pure
methane above 288 K and natural-gas mixtures are NOT validated (no measured
data were accessible).  CO2 is checked at Q1 only (−1.2 %); without the
dissolved-gas correction (before v1.19.1) CO2 formation pressures came out
15–30 % low, about 2 K.

Scope and limitations
---------------------
- Liquid water–hydrate–gas (Lw-H-V and Lw-H-Lhc) only, T ≥ 273.15 K; the
  ice branch is not implemented.
- Pure water, no inhibitors.  The hydrate formers dissolved in the water
  lower its activity (x_w < 1, Henry's law with NIST constants, no Poynting
  correction); this matters most for CO2 (x_CO2 ≈ 0.02–0.03).
- Structure H and heavy guests (C5+) are not hydrate formers here; they only
  dilute the gas and lower the formers' fugacities.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from functools import lru_cache

_R = 8.314462618
_T0 = 273.15
_ATM = 101_325.0
_H_FUS_ICE = 6009.5           # J/mol
_DV_FUS_ICE = 1.6e-6          # m³/mol, v_liquid − v_ice (liquid is denser)


@dataclass(frozen=True)
class _Structure:
    name: str
    nu_small: float
    nu_large: float
    dmu0: float                # J/mol
    dh0_ice: float             # J/mol
    dv0_ice: float             # m³/mol


STRUCTURES = (
    _Structure("sI", 1.0 / 23.0, 3.0 / 23.0, 1264.0, 1151.0, 3.0e-6),
    _Structure("sII", 2.0 / 17.0, 1.0 / 17.0, 883.0, 808.0, 3.4e-6),
)

# Munck et al. (1988): (A×10³ K/atm, B K) for (small, large) cavities.
_MUNCK = {
    # CAS:        sI small        sI large         sII small        sII large
    "74-82-8":   ((0.7228, 3187), (23.35, 2653), (0.2207, 3453), (100.0, 1916)),   # methane
    "74-84-0":   ((0.0, 0.0),     (3.039, 3861), (0.0, 0.0),     (240.0, 2967)),   # ethane
    "74-98-6":   ((0.0, 0.0),     (0.0, 0.0),    (0.0, 0.0),     (5.455, 4638)),   # propane
    "75-28-5":   ((0.0, 0.0),     (0.0, 0.0),    (0.0, 0.0),     (189.3, 3800)),   # isobutane
    "106-97-8":  ((0.0, 0.0),     (0.0, 0.0),    (0.0, 0.0),     (30.51, 3699)),   # n-butane
    "7727-37-9": ((1.617, 2905),  (6.078, 2431), (0.1742, 3082), (18.00, 1728)),   # nitrogen
    "124-38-9":  ((0.2474, 3410), (42.46, 2813), (0.0845, 3615), (851.0, 2025)),   # CO2
}
HYDRATE_FORMER_CAS = frozenset(_MUNCK)

# Henry's-law solubility in water, NIST Chemistry WebBook (R. Sander's
# compilation, NIST SRD 69): k_H(T) = k°_H·exp(D·(1/T − 1/298.15)),
# k°_H in mol/(kg·bar), D = d ln k_H / d(1/T) in K.  For each gas the most
# frequent entry that has a temperature dependence.
_HENRY = {
    "74-82-8": (0.0014, 1600.0),     # methane
    "74-84-0": (0.0019, 2300.0),     # ethane
    "74-98-6": (0.0015, 2700.0),     # propane
    "75-28-5": (0.00081, 2700.0),    # isobutane
    "106-97-8": (0.0012, 3100.0),    # n-butane
    "7727-37-9": (0.00065, 1300.0),  # nitrogen
    "124-38-9": (0.034, 2400.0),     # CO2
}
_M_WATER = 0.018015                   # kg/mol


def water_mole_fraction(T: float, fug_atm: dict[str, float]) -> float:
    """Mole fraction of water in the aqueous liquid with the dissolved hydrate
    formers: molality m_j = k_H,j(T)·f_j (Henry), x_w = (1/M_w) / (1/M_w + Σ m_j).
    The Poynting factor exp(−v̄_j^∞·P/RT) is neglected (no sourced partial
    molar volumes); it would lower the solubility by about 5–7 % at 5 MPa."""
    m_total = 0.0
    for cas, f_atm in fug_atm.items():
        k0, D = _HENRY.get(cas, (0.0, 0.0))
        m_total += k0 * math.exp(D * (1.0 / T - 1.0 / 298.15)) * f_atm * 1.01325
    return (1.0 / _M_WATER) / (1.0 / _M_WATER + m_total)


def _langmuir(cas: str, structure: int, cavity: int, T: float) -> float:
    """Langmuir constant (1/atm) of a guest in a cavity (0 small, 1 large)."""
    A, B = _MUNCK[cas][2 * structure + cavity]
    if A == 0.0:
        return 0.0
    return A * 1e-3 / T * math.exp(B / T)


def _dmu_water_liquid(st: _Structure, T: float, P: float, x_w: float = 1.0) -> float:
    """Δμ_w^{β−L}/RT on the liquid-water basis."""
    dh0_l = st.dh0_ice - _H_FUS_ICE
    dv_l = st.dv0_ice + _DV_FUS_ICE

    def dh(Tp: float) -> float:
        d = Tp - _T0
        return dh0_l - 38.12 * d + 0.5 * 0.141 * d * d

    # ∫ Δh/(R T²) dT by Simpson's rule (smooth integrand)
    n = 40
    a, b = _T0, T
    hstep = (b - a) / n
    s = dh(a) / a**2 + dh(b) / b**2
    for i in range(1, n):
        Tp = a + i * hstep
        s += (4 if i % 2 else 2) * dh(Tp) / Tp**2
    integral = s * hstep / 3.0 / _R
    return st.dmu0 / (_R * _T0) - integral + dv_l * P / (_R * T) - math.log(x_w)


def _dmu_water_hydrate(st_index: int, T: float, fug_atm: dict[str, float]) -> float:
    """Δμ_w^{β−H}/RT from cavity occupancies."""
    st = STRUCTURES[st_index]
    total = 0.0
    for cavity, nu in ((0, st.nu_small), (1, st.nu_large)):
        cf = sum(_langmuir(cas, st_index, cavity, T) * f for cas, f in fug_atm.items())
        total += nu * math.log1p(cf)          # −ν ln(1 − Σθ) = ν ln(1 + Σ C f)
    return total


def _former_fugacities(component_names, zs, T: float, P: float, eos_name: str) -> dict[str, float]:
    """Fugacities (atm) of the hydrate formers in the dry hydrocarbon fluid."""
    from .compositional_fluid import _get_flash_obj
    flash_obj = _get_flash_obj(tuple(component_names), eos_name)
    consts = flash_obj.constants
    try:
        res = flash_obj.flash(T=T, P=P, zs=list(zs))
        phase = res.gas if res.gas is not None else res.phases[0]
        fugs = phase.fugacities()
    except ZeroDivisionError:                 # single component: no stability test
        phase = flash_obj.gas.to(T=T, P=P, zs=list(zs))
        phase.constants, phase.correlations = flash_obj.gas.constants, flash_obj.gas.correlations
        fugs = phase.fugacities()
    return {cas: f / _ATM for cas, f in zip(consts.CASs, fugs) if cas in HYDRATE_FORMER_CAS and f > 0.0}


@dataclass(frozen=True)
class HydratePoint:
    temperature_c: float
    pressure_pa: float
    structure: str


def _dry_basis(component_names, zs):
    from .free_water import find_water_index
    w = find_water_index(tuple(component_names))
    names = tuple(n for i, n in enumerate(component_names) if i != w)
    z = [x for i, x in enumerate(zs) if i != w]
    s = sum(z)
    return names, tuple(x / s for x in z)


def _solve_formation_pressure(names: tuple[str, ...], zs: tuple[float, ...], T: float, eos_name: str,
                              p_guess: float | None = None):
    """(P_eq, structure) with the lowest formation pressure, or None.

    With ``p_guess`` (e.g. the previous point of a curve) the root is
    bracketed by ×1.25 steps around it; otherwise by stepping up from 0.2 bar
    (at very high pressure the Δv·P term can make the empty lattice stable
    again, so the ends of the range may share a sign).  Regula falsi
    (Illinois) refines the bracket to 1e-5 in ln P.
    """
    from .compositional_fluid import _get_flash_obj
    cas = _get_flash_obj(names, eos_name).constants.CASs
    if not any(c in HYDRATE_FORMER_CAS and z > 0.0 for c, z in zip(cas, zs)):
        return None
    ln_lo_lim, ln_hi_lim = math.log(2e4), math.log(1.5e8)
    step = math.log(1.25)
    best = None
    sii_formers = {"74-98-6", "75-28-5", "106-97-8"}          # propane, butanes: sII only
    order = (1, 0) if any(c in sii_formers and z > 0.0 for c, z in zip(cas, zs)) else (0, 1)
    for si in order:
        st = STRUCTURES[si]

        def F(lnP: float, st=st, si=si) -> float:
            P = math.exp(lnP)
            fug = _former_fugacities(names, zs, T, P, eos_name)
            return _dmu_water_liquid(st, T, P, water_mole_fraction(T, fug)) - _dmu_water_hydrate(si, T, fug)

        if best is not None and F(math.log(best[0])) > 0.0:
            continue                          # not stable at the other structure's pressure
        bracket = None
        guess = best[0] if best is not None else p_guess
        if guess is not None:
            a = min(max(math.log(guess), ln_lo_lim), ln_hi_lim)
            f_a = F(a)
            if f_a > 0.0:                      # below the root: step up
                while a < ln_hi_lim and bracket is None:
                    b = a + step
                    f_b = F(b)
                    if f_b <= 0.0:
                        bracket = (a, b, f_a, f_b)
                    a, f_a = b, f_b
            else:                              # above the root: step down
                b, f_b = a, f_a
                while b > ln_lo_lim and bracket is None:
                    a = b - step
                    f_a = F(a)
                    if f_a > 0.0:
                        bracket = (a, b, f_a, f_b)
                    b, f_b = a, f_a
        if bracket is None:
            a = ln_lo_lim
            f_a = F(a)
            while a < ln_hi_lim:
                b = a + step
                f_b = F(b)
                if f_a > 0.0 >= f_b:
                    bracket = (a, b, f_a, f_b)
                    break
                a, f_a = b, f_b
        if bracket is None:
            continue                          # this structure does not form in range
        lo, hi, f_lo, f_hi = bracket
        side = 0
        x = 0.5 * (lo + hi)
        for _ in range(60):
            x = hi - f_hi * (hi - lo) / (f_hi - f_lo) if f_hi != f_lo else 0.5 * (lo + hi)
            fx = F(x)
            if abs(fx) < 1e-12 or hi - lo < 1e-5:
                break
            if fx > 0.0:
                lo, f_lo = x, fx
                if side == 1:
                    f_hi *= 0.5
                side = 1
            else:
                hi, f_hi = x, fx
                if side == -1:
                    f_lo *= 0.5
                side = -1
        P_eq = math.exp(x)
        if best is None or P_eq < best[0]:
            best = (P_eq, st.name)
    return best


@lru_cache(maxsize=4096)
def _formation_pressure(names: tuple[str, ...], zs: tuple[float, ...], T: float, eos_name: str):
    return _solve_formation_pressure(names, zs, T, eos_name)


def hydrate_formation_pressure(component_names, zs, temperature_c: float, eos_name: str = "PR") -> HydratePoint | None:
    """Pressure (Pa) above which hydrate forms at ``temperature_c`` (≥ 0 °C),
    and the stable structure; None if no former is present or it does not
    form below 1500 bar."""
    names, z = _dry_basis(component_names, zs)
    T = temperature_c + 273.15
    if T < _T0 - 1e-9:
        raise ValueError("The ice branch (T < 0 °C) is not implemented.")
    z_r = tuple(round(x, 8) for x in z)
    res = _formation_pressure(names, z_r, round(T, 4), eos_name.upper())
    if res is None:
        return None
    return HydratePoint(temperature_c, res[0], res[1])


@lru_cache(maxsize=4096)
def _formation_temperature(names: tuple[str, ...], zs: tuple[float, ...], P: float, eos_name: str):
    """Highest T (K) at which some structure is stable at pressure P, or None."""
    from .compositional_fluid import _get_flash_obj
    cas = _get_flash_obj(names, eos_name).constants.CASs
    if not any(c in HYDRATE_FORMER_CAS and z > 0.0 for c, z in zip(cas, zs)):
        return None
    best = None
    T_lo0, T_hi0 = _T0, _T0 + 45.0
    for si, st in enumerate(STRUCTURES):
        def G(T: float) -> float:
            fug = _former_fugacities(names, zs, T, P, eos_name)
            return _dmu_water_liquid(st, T, P, water_mole_fraction(T, fug)) - _dmu_water_hydrate(si, T, fug)

        g_lo = G(T_lo0)
        if g_lo > 0.0:
            continue                          # not stable even at 0 °C
        g_hi = G(T_hi0)
        if g_hi <= 0.0:
            T_eq = T_hi0                      # still stable at 45 °C: report the cap
        else:
            # Regula falsi (Illinois): G is smooth in T, a few steps suffice.
            lo, hi, g_l, g_h = T_lo0, T_hi0, g_lo, g_hi
            T_eq = 0.5 * (lo + hi)
            side = 0
            for _ in range(60):
                T_eq = hi - g_h * (hi - lo) / (g_h - g_l)
                g = G(T_eq)
                if abs(g) < 1e-9 or hi - lo < 1e-4:
                    break
                if g <= 0.0:
                    lo, g_l = T_eq, g
                    if side == -1:
                        g_h *= 0.5
                    side = -1
                else:
                    hi, g_h = T_eq, g
                    if side == 1:
                        g_l *= 0.5
                    side = 1
        if best is None or T_eq > best[0]:
            best = (T_eq, st.name)
    return best


def hydrate_formation_temperature(component_names, zs, pressure_pa: float, eos_name: str = "PR") -> HydratePoint | None:
    """Temperature (°C) below which hydrate forms at ``pressure_pa``, and the
    stable structure; None when no hydrate is stable above 0 °C at that
    pressure (45 °C is reported as a cap)."""
    names, z = _dry_basis(component_names, zs)
    z_r = tuple(round(x, 8) for x in z)
    res = _formation_temperature(names, z_r, float(round(pressure_pa, -2)), eos_name.upper())
    if res is None:
        return None
    return HydratePoint(res[0] - 273.15, pressure_pa, res[1])


def hydrate_curve(component_names, zs, t_min_c: float = 0.0, t_max_c: float = 30.0,
                  n: int = 31, eos_name: str = "PR", p_max_pa: float = 1.0e8) -> list[HydratePoint]:
    """Formation curve: points (T, P_formation, structure) from t_min to t_max,
    stopping after the first point above ``p_max_pa`` (the curve is steep, and
    points far above a plot's range only cost time).  Each point starts its
    search from the previous one (warm start)."""
    names, z = _dry_basis(component_names, zs)
    z_r = tuple(round(x, 8) for x in z)
    pts: list[HydratePoint] = []
    guess = None
    for i in range(n):
        T = t_min_c + (t_max_c - t_min_c) * i / max(n - 1, 1)
        res = _solve_formation_pressure(names, z_r, max(T, 0.0) + 273.15, eos_name.upper(), guess)
        if res is None:
            continue
        guess = res[0]
        pts.append(HydratePoint(T, res[0], res[1]))
        if res[0] > p_max_pa:
            break
    return pts


# ── Gas-gravity method (Katz chart) for black-oil fluids ──────────────────────
#
# A black-oil fluid knows only its gas gravity, not a composition, so its
# hydrate check uses the classic gas-gravity chart of Katz (1945), reproduced
# in the GPSA Engineering Data Book.  Chart readings (pressure in psia at
# 30–80 °F for gravities 0.6–1.0) as digitised in F. Janna, "Correlation for
# predicting hydrate formation", B.Eng. dissertation, Universiti Teknologi
# PETRONAS (2015), Tables 4.1a–4.5a ("Experiment" column = chart values).
# Interpolation is linear in ln P and in gravity.  Below 0.6 the chart is
# anchored at pure methane (0.5537) by the van der Waals–Platteeuw curve
# above.  The chart is an engineering estimate: it does not see the propane
# or CO2 content that a composition-based prediction would.
_KATZ_T_F = (30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0)
_KATZ_P_PSIA = {
    0.6: (130, 185, 258, 359, 510, 709, 990, 1400, 2300, 4220, 7800),
    0.7: (85, 123, 175, 249, 358, 510, 729, 1063, 1937, 3539, 6500),
    0.8: (68, 98, 143, 209, 299, 423, 611, 940, 1709, 3080, 5650),
    0.9: (55, 77, 115, 170, 240, 350, 510, 840, 1500, 2700, 5000),
    1.0: (45, 62, 90, 127.13, 190, 299, 424.53, 720, 1280.25, 2400, 4265),
}
_METHANE_GRAVITY = 16.043 / 28.964
_PSI = 6894.757


def _katz_lnP_curve(gas_gravity: float) -> list[float]:
    """ln P (Pa) at the chart temperatures for a gas gravity (interpolated)."""
    grav = sorted(_KATZ_P_PSIA)
    g = min(max(gas_gravity, _METHANE_GRAVITY), grav[-1])
    if g >= grav[0]:
        for g0, g1 in zip(grav, grav[1:]):
            if g0 <= g <= g1:
                w = (g - g0) / (g1 - g0)
                return [(1 - w) * math.log(a * _PSI) + w * math.log(b * _PSI)
                        for a, b in zip(_KATZ_P_PSIA[g0], _KATZ_P_PSIA[g1])]
    # between pure methane (vdW-P) and the 0.6 chart line
    w = (g - _METHANE_GRAVITY) / (grav[0] - _METHANE_GRAVITY)
    out = []
    for t_f, p06 in zip(_KATZ_T_F, _KATZ_P_PSIA[grav[0]]):
        t_c = max((t_f - 32.0) / 1.8, 0.0)
        pm = hydrate_formation_pressure(["methane"], [1.0], t_c)
        ln_m = math.log(pm.pressure_pa) if pm is not None else math.log(p06 * _PSI)
        out.append((1 - w) * ln_m + w * math.log(p06 * _PSI))
    return out


def katz_hydrate_pressure(gas_gravity: float, temperature_c: float) -> float:
    """Hydrate formation pressure (Pa) from the gas-gravity chart."""
    t_f = temperature_c * 1.8 + 32.0
    lnp = _katz_lnP_curve(gas_gravity)
    Ts = _KATZ_T_F
    i = 0 if t_f <= Ts[0] else (len(Ts) - 2 if t_f >= Ts[-1] else
                                 max(k for k in range(len(Ts) - 1) if Ts[k] <= t_f))
    w = (t_f - Ts[i]) / (Ts[i + 1] - Ts[i])
    return math.exp(lnp[i] + w * (lnp[i + 1] - lnp[i]))


def katz_hydrate_temperature(gas_gravity: float, pressure_pa: float) -> HydratePoint | None:
    """Hydrate formation temperature (°C) at ``pressure_pa`` from the chart;
    None below the 0 °C pressure, capped at the chart's 80 °F (26.7 °C)."""
    lnp = _katz_lnP_curve(gas_gravity)
    ln_target = math.log(max(pressure_pa, 1.0))
    Ts = _KATZ_T_F
    p0 = katz_hydrate_pressure(gas_gravity, 0.0)
    if pressure_pa < p0:
        return None
    if ln_target >= lnp[-1]:
        return HydratePoint((Ts[-1] - 32.0) / 1.8, pressure_pa, "gas-gravity chart (cap)")
    for k in range(len(Ts) - 1):
        if lnp[k] <= ln_target <= lnp[k + 1]:
            t_f = Ts[k] + (ln_target - lnp[k]) / (lnp[k + 1] - lnp[k]) * (Ts[k + 1] - Ts[k])
            return HydratePoint(max((t_f - 32.0) / 1.8, 0.0), pressure_pa, "gas-gravity chart")
    return HydratePoint(0.0, pressure_pa, "gas-gravity chart")
