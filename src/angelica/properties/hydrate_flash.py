"""
Hydrate as a phase: how much hydrate forms (roadmap item 2b.2).

At fixed T and P, a feed of water and hydrocarbons is split into a fluid
(hydrocarbon vapour/liquid, free water, water vapour) and a hydrate phase
of one structure.  The equilibrium is the minimum of the Gibbs energy; with
the van der Waals–Platteeuw model it is reached where the chemical
potential of water is the same in the hydrate and in the fluid:

    G(ξ) = Δμ_w^{β−L}(a_w) − Δμ_w^{β−H}(f) = 0

ξ is the water in the hydrate (mol per mol feed).  The hydrate phase has one
composition, set by the cavity occupancies θ_mj at the final fugacities
f_j of the fluid: guest j per water molecule g_j = Σ_m ν_m·θ_mj.  The fluid
keeps z_j − ξ·g_j of each former and z_w − ξ of the water, so the mass
balances close by construction.  As hydrate forms, three things can stop
it:

  - ``"equilibrium"``: the formers' fugacities fall (the gas loses its best
    formers) until G = 0, with free water left;
  - ``"water"``: the free water is used up; the gas's water then sets the
    water activity y_w·P/P_sat(T) (hydrate–vapour equilibrium, roadmap
    2b.1) and G = 0 at an undersaturated water content.  If the
    hydrocarbon has no vapour (an oil, or a dense gas the EOS calls liquid)
    the free-water model gives it no dissolved water: all the water goes
    into the hydrate and G stays < 0;
  - ``"formers"``: a former runs out with G still < 0 (e.g. pure methane,
    whose fugacity does not fall as it is consumed).

``"none"``: no hydrate is stable at the feed (G(0) ≥ 0).  The structure is
the one most stable at the feed (most negative G(0)); if the other is also
stable at the end, ``other_structure_stable`` is set (a two-hydrate
equilibrium is not computed, so the result is then not the Gibbs minimum;
near where the other structure takes over the amount varies steeply).

Water in the fluid is distributed with Angelica's free-water model
(immiscible water, Raoult for the vapour).  Its limit shows here: a dense
gas that the EOS labels liquid above some pressure carries no water vapour,
and the amount of hydrate jumps there (C1/C2/C3 with 0.05 % water at
15 °C: 3.79e-4 mol at 130 bar, 5.0e-4 at 140 bar).  Scope: T ≥ 0 °C (no ice phase,
as in the vdW-P curve), no inhibitors.  Black-oil fluids have no
composition: the amount of hydrate is a compositional calculation only.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

from .hydrate import (
    _T0,
    HYDRATE_FORMER_CAS,
    STRUCTURES,
    _dmu_water_hydrate,
    _dmu_water_liquid,
    _former_fugacities,
    _langmuir,
    water_mole_fraction,
)

_EXHAUSTED = 1e-12           # a former's amount (mol/mol feed) below which it is used up
_MAX_NEWTON = 60            # Newton iterations per ξ


@dataclass(frozen=True)
class HydratePartition:
    """Split of 1 mol of feed at (T, P) into fluid and hydrate."""
    temperature_c: float
    pressure_pa: float
    limited_by: str                        # "none", "equilibrium", "water", "formers"
    structure: str | None
    hydrate_water_mol: float               # water in the hydrate, mol per mol feed
    hydrate_guests_mol: dict[str, float]   # name → mol per mol feed
    fluid_mol: dict[str, float]            # name → mol per mol feed left in the fluid (water included)
    free_water_mol: float                  # liquid water in the fluid, mol per mol feed
    vapour_water_mol: float                # water vapour in the fluid, mol per mol feed
    driving_force: float                   # G/RT at the result (0 at equilibrium, < 0 if a former ran out)
    other_structure_stable: bool = False
    hydration_number: float | None = None  # water molecules per guest in the hydrate
    notes: tuple[str, ...] = field(default_factory=tuple)

    @property
    def hydrate_mol(self) -> float:
        return self.hydrate_water_mol + sum(self.hydrate_guests_mol.values())


def _guest_ratio(si: int, T: float, fug: dict[str, float]) -> dict[str, float]:
    """Guests per water molecule in structure ``si`` at fugacities ``fug`` (atm)."""
    st = STRUCTURES[si]
    out: dict[str, float] = {}
    for cavity, nu in ((0, st.nu_small), (1, st.nu_large)):
        cf = {cas: _langmuir(cas, si, cavity, T) * f for cas, f in fug.items()}
        denom = 1.0 + sum(cf.values())
        for cas, c in cf.items():
            if c > 0.0:
                out[cas] = out.get(cas, 0.0) + nu * c / denom
    return out


def _solve(A: list[list[float]], b: list[float]) -> list[float] | None:
    """Gaussian elimination with partial pivoting (small systems)."""
    n = len(b)
    M = [row[:] + [bi] for row, bi in zip(A, b)]
    for c in range(n):
        p = max(range(c, n), key=lambda i: abs(M[i][c]))
        if abs(M[p][c]) < 1e-300:
            return None
        M[c], M[p] = M[p], M[c]
        for i in range(c + 1, n):
            f = M[i][c] / M[c][c]
            for k in range(c, n + 1):
                M[i][k] -= f * M[c][k]
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        x[i] = (M[i][n] - sum(M[i][k] * x[k] for k in range(i + 1, n))) / M[i][i]
    return x


def hydrate_partition(component_names, zs, temperature_c: float, pressure_pa: float,
                      eos_name: str = "PR") -> HydratePartition:
    """Equilibrium amount of hydrate formed by 1 mol of feed ``zs`` (with a
    water component) at ``temperature_c`` and ``pressure_pa``."""
    from .compositional_fluid import _get_flash_obj, _hc_phase_compositions
    from .free_water import find_water_index, free_water_split, water_psat_pa

    names = tuple(component_names)
    w = find_water_index(names)
    if w is None:
        raise ValueError("hydrate_partition needs a water component in the feed")
    if not (math.isfinite(temperature_c) and math.isfinite(pressure_pa) and pressure_pa > 0.0):
        raise ValueError(f"need a finite temperature and a pressure > 0 (got {temperature_c} °C, {pressure_pa} Pa)")
    T = temperature_c + 273.15
    if T < _T0 - 1e-9:
        raise ValueError("The ice branch (T < 0 °C) is not implemented.")
    z = [float(x) for x in zs]
    if len(z) != len(names) or any(not math.isfinite(x) or x < 0.0 for x in z) or sum(z) <= 0.0:
        raise ValueError("zs must be finite, ≥ 0 and as long as component_names")
    s = sum(z)
    z = [x / s for x in z]
    eos = eos_name.upper()
    hc_names = tuple(n for i, n in enumerate(names) if i != w)
    cas = _get_flash_obj(hc_names, eos).constants.CASs
    name_of = dict(zip(cas, hc_names))
    z_hc = [x for i, x in enumerate(z) if i != w]
    z_w = z[w]
    P = float(pressure_pa)
    psat = water_psat_pa(T)

    def fluid_state(n_hc: list[float], n_w: float):
        """Fugacities, water activity, (free, vapour) water of the fluid."""
        tot = sum(n_hc)
        zz = tuple(x / tot for x in n_hc)
        fug = _former_fugacities(hc_names, zz, T, P, eos)
        vf = _hc_phase_compositions(hc_names, P, temperature_c, tuple(round(x, 14) for x in zz), eos)[2]
        n_wl, n_wv = free_water_split(n_w, vf, tot, T, P)
        v_hc = vf * tot
        if n_wl > 1e-15 or v_hc <= 0.0:
            a_w = water_mole_fraction(T, fug)
        else:
            y_w = n_wv / (n_wv + v_hc)
            a_w = min(water_mole_fraction(T, fug), y_w * P / psat)
        return fug, a_w, n_wl, n_wv

    def g_value(si: int, fug, a_w: float) -> float:
        if a_w <= 0.0:
            return math.inf
        return _dmu_water_liquid(STRUCTURES[si], T, P, a_w) - _dmu_water_hydrate(si, T, fug)

    fug0, a0, wl0, wv0 = fluid_state(z_hc, z_w)
    g0 = [g_value(si, fug0, a0) for si in range(len(STRUCTURES))]
    no_hydrate = HydratePartition(temperature_c, P, "none", None, 0.0, {},
                                  dict(zip(names, z)), wl0, wv0, min(g0))
    if not any(c in HYDRATE_FORMER_CAS and x > 0.0 for c, x in zip(cas, z_hc)) or z_w <= 0.0 or min(g0) >= 0.0:
        return no_hydrate
    si = min(range(len(STRUCTURES)), key=lambda i: g0[i])

    warm: dict[str, list[float]] = {}

    formers = [j for j, (c, x) in enumerate(zip(cas, z_hc))
               if x > 0.0 and c in HYDRATE_FORMER_CAS
               and any(_langmuir(c, si, m, T) > 0.0 for m in (0, 1))]

    def evaluate(xi: float, u: list[float]):
        n = list(z_hc)
        for k, j in enumerate(formers):
            n[j] = math.exp(u[k])
        fug, a_w, n_wl, n_wv = fluid_state(n, z_w - xi)
        ratio = _guest_ratio(si, T, fug)
        r = [(n[j] + xi * ratio.get(cas[j], 0.0) - z_hc[j]) / z_hc[j] for j in formers]
        return r, (n, ratio, fug, a_w, n_wl, n_wv)

    def state_at(xi: float):
        """Fluid at hydrate water ξ: solve n_j + ξ·g_j(f(n)) = z_j for the
        formers by Newton in ln n_j (finite-difference Jacobian, step limited,
        backtracking).  None when a former runs out (no solution with n > 0:
        n_j → 0 while z_j < ξ·g_j)."""
        u = list(warm.get("u", [math.log(z_hc[j]) for j in formers]))
        r, st = evaluate(xi, u)
        for _ in range(_MAX_NEWTON):
            norm = max(abs(x) for x in r)
            if norm < 1e-13:
                break
            jac = []
            for k in range(len(u)):
                up = list(u)
                up[k] += 1e-6
                rp, _ = evaluate(xi, up)
                jac.append([(a - b) / 1e-6 for a, b in zip(rp, r)])      # column k
            A = [[jac[k][i] for k in range(len(u))] for i in range(len(u))]
            step = _solve(A, [-x for x in r])
            if step is None:
                return None
            big = max(abs(x) for x in step)
            if big > 8.0:
                step = [x * 8.0 / big for x in step]
            lam = 1.0
            for _ in range(20):
                u_try = [a + lam * b for a, b in zip(u, step)]
                r_try, st_try = evaluate(xi, u_try)
                if max(abs(x) for x in r_try) < norm:
                    break
                lam *= 0.5
            u, r, st = u_try, r_try, st_try
            n_now, ratio_now = st[0], st[1]
            if any(n_now[j] < 1e-7 * z_hc[j] and xi * ratio_now.get(cas[j], 0.0) > z_hc[j] * (1.0 + 1e-7)
                   for j in formers) or any(ui - math.log(z_hc[j]) < -30.0 for ui, j in zip(u, formers)):
                return None                            # a former is used up: z_j < ξ·g_j as n_j → 0
        if max(abs(x) for x in r) > 1e-9:
            return None
        warm["u"] = u
        n, ratio, fug, a_w, n_wl, n_wv = st
        n = [zj - xi * ratio.get(c, 0.0) for zj, c in zip(z_hc, cas)]     # the balance, exactly
        return n, ratio, fug, a_w, n_wl, n_wv, g_value(si, fug, a_w)

    # Bracket ξ in [0, z_w): G < 0 (more hydrate) below the root, G ≥ 0 or a
    # former exhausted above it.  The bracket is first narrowed by decades
    # (a feed with little gas forms little hydrate), then bisected to 1e-10
    # of its upper end.
    lo, hi = 0.0, z_w * (1.0 - 1e-12)
    st_lo = state_at(lo)
    st_hi = state_at(hi)
    if st_lo is None:                                  # cannot happen: ξ = 0 is the feed
        return no_hydrate
    if st_hi is not None and st_hi[-1] < 0.0:
        # Still stable with all the water in the hydrate: this happens when
        # the hydrocarbon has no vapour (the free-water model's liquid
        # hydrocarbon holds no water, so a_w stays at the aqueous value).
        lo, st_lo = hi, st_hi
    else:
        x = hi
        for _ in range(40):
            x *= 0.1
            st = state_at(x)
            if st is not None and st[-1] < 0.0:
                lo, st_lo = x, st
                break
            hi, st_hi = x, st
        for _ in range(200):
            if hi - lo < 1e-10 * hi:
                break
            mid = 0.5 * (lo + hi)
            st = state_at(mid)
            if st is not None and st[-1] < 0.0:
                lo, st_lo = mid, st
            else:
                hi, st_hi = mid, st
    xi, st = lo, st_lo
    if st_hi is not None:
        # G changes sign in [lo, hi]: one secant step puts G at 0 (the root
        # of a smooth function; the bracket is 1e-10 of ξ wide).
        g_lo, g_hi = st_lo[-1], st_hi[-1]
        x = lo + (hi - lo) * (-g_lo) / (g_hi - g_lo) if g_hi != g_lo else lo
        st_x = state_at(x)
        if st_x is not None:
            xi, st = x, st_x
    n, ratio, fug, a_w, n_wl, n_wv, g = st
    guests_total = sum(ratio.values())
    # The bisection ends at G = 0 (equilibrium, or water used up) or at the
    # edge of the feasible ξ with G still < 0 (a former used up).
    if st_hi is None:
        limited = "formers"
    elif n_wl <= 1e-9 * z_w:
        limited = "water"
    else:
        limited = "equilibrium"
    other = [i for i in range(len(STRUCTURES)) if i != si and g_value(i, fug, a_w) < 0.0]
    fluid = {name_of[c]: x for c, x in zip(cas, n)}
    fluid[names[w]] = z_w - xi
    return HydratePartition(
        temperature_c=temperature_c, pressure_pa=P, limited_by=limited, structure=STRUCTURES[si].name,
        hydrate_water_mol=xi,
        hydrate_guests_mol={name_of[c]: xi * r for c, r in ratio.items() if xi * r > 0.0},
        fluid_mol=fluid, free_water_mol=n_wl, vapour_water_mol=n_wv, driving_force=g,
        other_structure_stable=bool(other),
        hydration_number=1.0 / guests_total if guests_total > 0.0 else None,
    )
