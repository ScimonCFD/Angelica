from __future__ import annotations

import math
from functools import cached_property, lru_cache

from .base import FluidModel
from .free_water import (
    WATER_MW,
    find_water_index,
    free_water_split,
    water_liquid_density_kg_m3,
    water_liquid_viscosity_pa_s,
)
from .surface_tension import water_surface_tension

# Flash-object cache keyed by (component_names, eos_name) — avoids recreating
# ChemicalConstantsPackage + CEOSGas/Liquid on every cache miss.
_FLASH_OBJS: dict = {}

# Supported EOS keys and their descriptions (shown in GUI combobox).
SUPPORTED_EOS: dict[str, str] = {
    "PR":     "Peng-Robinson (1976) — industry standard for gas/condensate",
    "PR78":   "Peng-Robinson (1978) — improved alpha for heavier HC",
    "PRSV":   "PR Stryjek-Vera — better for polar components (CO₂, H₂S, H₂O)",
    "SRK":    "Soave-Redlich-Kwong (1972)",
    "APISRK": "API-SRK — API Technical Data Book standard for petroleum",
}

_EOS_CLASS_MAP = {
    "PR":     "PRMIX",
    "PR78":   "PR78MIX",
    "PRSV":   "PRSVMIX",
    "SRK":    "SRKMIX",
    "APISRK": "APISRKMIX",
}


def _get_flash_obj(component_names: tuple[str, ...], eos_name: str):
    key = (component_names, eos_name)
    if key not in _FLASH_OBJS:
        import warnings

        import thermo.eos_mix as _eos_mix
        from thermo import ChemicalConstantsPackage
        from thermo.flash import FlashVL
        from thermo.phases import CEOSGas, CEOSLiquid

        from .phase_envelope import _build_kij_matrix
        eos_cls = getattr(_eos_mix, _EOS_CLASS_MAP.get(eos_name, "PRMIX"))
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            constants, props = ChemicalConstantsPackage.from_IDs(list(component_names))
            kijs = _build_kij_matrix(constants)
            eos_kw = {"Tcs": constants.Tcs, "Pcs": constants.Pcs, "omegas": constants.omegas, "kijs": kijs}
            gas_phase = CEOSGas(eos_cls, eos_kw, HeatCapacityGases=props.HeatCapacityGases)
            liq_phase = CEOSLiquid(eos_cls, eos_kw, HeatCapacityGases=props.HeatCapacityGases)
            flash_obj = FlashVL(constants, props, liquid=liq_phase, gas=gas_phase)
        _FLASH_OBJS[key] = flash_obj
    return _FLASH_OBJS[key]


# Flash results are memoised on rounded (P, T, z).  The rounding must be far
# finer than the outer-loop tolerances: with T rounded to 0.1 K a two-phase
# density moves by ~5e-4 per rounding step (5x the default density tolerance),
# and the outer loop can lock into a two-state limit cycle.
_P_STEP_PA = 10.0
_T_DIGITS = 3
_ZS_DIGITS = 6


def _cache_key_pt(P: float, T: float) -> tuple[float, float]:
    return float(round(P / _P_STEP_PA) * _P_STEP_PA), round(T, _T_DIGITS)


@lru_cache(maxsize=8192)
def _tp_flash(
    component_names: tuple[str, ...],
    pressure_pa: float,
    temperature_c: float,
    zs: tuple[float, ...],
    eos_name: str = "PR",
):
    """thermo PT-flash result, shared by every property of a state (density,
    phase split, enthalpy, phase compositions) so each state is flashed once.
    Exceptions (single-component stability failures) are not cached and
    propagate to the caller's fallback."""
    flash_obj = _get_flash_obj(component_names, eos_name)
    return flash_obj.flash(T=temperature_c + 273.15, P=pressure_pa, zs=list(zs))


# ── Viscosity ─────────────────────────────────────────────────────────────────
# thermo evaluates gas viscosity with low-pressure mixing rules only (no
# pressure dependence), returns a 1e-5 Pa·s placeholder when a dense
# supercritical phase is labelled "liquid", and its McAdams two-phase blend
# takes the mass fraction of phase 0 as the gas quality even when phase 0 is
# the liquid.  Viscosity is therefore assembled here:
#   gas-like phase : thermo low-pressure mixture value + Jossi-Stiel-Thodos
#                    dense-gas residual on the EOS reduced density
#   liquid phase   : thermo liquid mixture value (DIPPR-based, accurate)
#   two-phase      : McAdams, 1/μ = x_g/μ_g + (1 − x_g)/μ_l, x_g = gas mass fraction

_ATM_PA = 101_325.0


def _jst_residual_pa_s(rho_r: float, Tpc: float, Ppc: float, MW: float) -> float:
    """Jossi-Stiel-Thodos (1962) dense-fluid residual viscosity μ − μ⁰ [Pa·s].

    ρ_r = V_pc / V (Kay pseudo-critical volume over EOS molar volume).  This is
    the pressure term of the Lohrenz-Bray-Clark (1964) correlation, valid for
    gases and dense supercritical fluids up to ρ_r ≈ 3.
    """
    poly = 0.1023 + rho_r * (0.023364 + rho_r * (0.058533 + rho_r * (-0.040758 + rho_r * 0.0093724)))
    xi = Tpc ** (1.0 / 6.0) * MW ** -0.5 * (Ppc / _ATM_PA) ** (-2.0 / 3.0)
    return max(poly ** 4 - 1.0e-4, 0.0) / xi * 1.0e-3   # cP → Pa·s


def _gas_like(phase, T_K: float, constants) -> bool:
    """Gas-labelled phases, and any phase above its pseudo-critical temperature."""
    from thermo.phases import CEOSGas
    if isinstance(phase, CEOSGas):
        return True
    Tpc = sum(z * tc for z, tc in zip(phase.zs, constants.Tcs))
    return T_K >= Tpc


def _phase_viscosity(phase, T_K: float, P_pa: float, constants, correlations,
                     gas_like: bool | None = None) -> float:
    """Viscosity of one EOS phase [Pa·s] (see module notes above)."""
    zs = list(phase.zs)
    if gas_like is None:
        gas_like = _gas_like(phase, T_K, constants)
    if not gas_like:
        return float(phase.mu())
    ws = phase.ws() if hasattr(phase, "ws") else None
    mu0 = float(correlations.ViscosityGasMixture.mixture_property(T_K, P_pa, zs, ws))
    Tpc = sum(z * v for z, v in zip(zs, constants.Tcs))
    Ppc = sum(z * v for z, v in zip(zs, constants.Pcs))
    Vpc = sum(z * v for z, v in zip(zs, constants.Vcs))
    MW = sum(z * v for z, v in zip(zs, constants.MWs))
    rho_r = Vpc / float(phase.V())
    return mu0 + _jst_residual_pa_s(rho_r, Tpc, Ppc, MW)


def _result_viscosity(res, T_K: float, P_pa: float, constants, correlations) -> float:
    """Bulk viscosity of a flash result: single phase, or McAdams blend."""
    phases = list(res.phases)
    if len(phases) == 1:
        return _phase_viscosity(phases[0], T_K, P_pa, constants, correlations)
    gas_like = [_gas_like(ph, T_K, constants) for ph in phases]
    if all(gas_like) or not any(gas_like):
        # thermo can label both phases alike near the critical point:
        # treat the lighter one as the vapour.
        i_gas = min(range(len(phases)), key=lambda i: phases[i].rho_mass())
        gas_like = [i == i_gas for i in range(len(phases))]
    mus = [_phase_viscosity(ph, T_K, P_pa, constants, correlations, g)
           for ph, g in zip(phases, gas_like)]
    betas_mass = list(res.betas_mass)
    inv = sum(b / mu for b, mu in zip(betas_mass, mus))
    return 1.0 / inv if inv > 0 else float(res.mu())


def _wilson_psat(Tc: float, Pc: float, omega: float, T_K: float) -> float:
    """Wilson K-value estimate for saturation pressure (Pa) of a single component."""
    if T_K >= Tc:
        return Pc
    return Pc * math.exp(5.373 * (1.0 + omega) * (1.0 - Tc / T_K))


_SIGMA_PLACEHOLDER = 0.020

# Parachors, (dyn/cm)^¼·cm³/mol, keyed by CAS number (Weinaug & Katz 1943;
# values as tabulated in Whitson & Brulé, *Phase Behavior*, SPE Monograph 20,
# 2000, and Reid, Prausnitz & Poling).  Components not listed use the Fanchi
# (1985) correlation P = −11.4 + 3.23·M − 0.0022·M².
_PARACHORS: dict[str, float] = {
    "74-82-8": 77.0,     # methane
    "74-84-0": 108.0,    # ethane
    "74-98-6": 150.3,    # propane
    "75-28-5": 181.5,    # isobutane
    "106-97-8": 189.9,   # n-butane
    "78-78-4": 225.0,    # isopentane
    "109-66-0": 231.5,   # n-pentane
    "110-54-3": 271.0,   # n-hexane
    "142-82-5": 312.5,   # n-heptane
    "111-65-9": 351.5,   # n-octane
    "111-84-2": 393.0,   # n-nonane
    "124-18-5": 433.5,   # n-decane
    "7727-37-9": 41.0,   # nitrogen
    "124-38-9": 78.0,    # carbon dioxide
    "7783-06-4": 80.1,   # hydrogen sulfide
}


def _parachors(constants) -> list[float]:
    out = []
    for cas, mw in zip(constants.CASs, constants.MWs):
        p = _PARACHORS.get(cas)
        out.append(p if p is not None else max(-11.4 + 3.23 * mw - 0.0022 * mw * mw, 10.0))
    return out


_R_GAS = 8.314462618


def _peneloux_molar_density(constants, phase) -> float:
    """Molar density (mol/cm³) with the Peneloux volume translation for PR,
    c_i = 0.50033·R·Tc/Pc·(0.25969 − Z_RA), Z_RA = 0.29056 − 0.08775·ω.

    Macleod-Sugden scales with density to the fourth power, and untranslated
    PR liquid densities are 3–8 % off for C5–C10 (σ off by 15–30 %); the shift
    brings pure n-alkane σ at 20 °C within 5 % of measured values.  It is used
    for σ only; the flow model keeps the EOS densities.
    """
    shift = sum(
        x * 0.50033 * _R_GAS * tc / pc * (0.25969 - (0.29056 - 0.08775 * w))
        for x, tc, pc, w in zip(phase.zs, constants.Tcs, constants.Pcs, constants.omegas)
    )
    V = float(phase.V()) - shift                                   # m³/mol
    return 1e-6 / V if V > 0.0 else float(phase.rho_mass()) * 1e-3 / float(phase.MW())


def _macleod_sugden_sigma(constants, liq, gas=None) -> float:
    """Gas–liquid surface tension (N/m) by Macleod-Sugden / Weinaug-Katz:

        σ^¼ = Σ P_i·(x_i·ρ_L/M_L − y_i·ρ_V/M_V)    (cgs: dyn/cm, mol/cm³)

    with Peneloux-translated molar densities.  The vapour term makes σ
    vanish at the critical point.  Without a gas phase (single-phase liquid)
    the vapour term is dropped.
    """
    par = _parachors(constants)
    s = sum(p * x for p, x in zip(par, liq.zs)) * _peneloux_molar_density(constants, liq)
    if gas is not None:
        s -= sum(p * y for p, y in zip(par, gas.zs)) * _peneloux_molar_density(constants, gas)
    if not math.isfinite(s) or s <= 0.0:
        return 1e-5
    return max(s ** 4 * 1e-3, 1e-5)                               # dyn/cm → N/m


@lru_cache(maxsize=2048)
def _hc_phase_split(
    component_names: tuple[str, ...],
    pressure_pa: float,
    temperature_c: float,
    zs: tuple[float, ...],
    eos_name: str = "PR",
) -> tuple[float, float, float, float, float, float, float, float]:
    """Per-phase properties of an HC-only PT flash.

    Returns ``(beta_g, rho_g, rho_l, MW_g, MW_l, mu_g, mu_l, sigma_l)``: the
    vapour mole fraction; for each phase, mass density (kg/m³), molar mass
    (g/mol) and viscosity (Pa·s); and the gas–liquid surface tension (N/m,
    Macleod-Sugden with parachors; 0.02 N/m placeholder when there is no
    liquid).  When two phases coexist the lighter one is the vapour
    (thermo can label both alike near the critical point).  A single-phase
    result repeats its properties in both slots, with ``beta_g`` 1 or 0.
    """
    flash_obj = _get_flash_obj(component_names, eos_name)
    const, corr = flash_obj.constants, flash_obj.correlations
    T_K = temperature_c + 273.15
    try:
        res = _tp_flash(component_names, pressure_pa, temperature_c, zs, eos_name)
    except ZeroDivisionError:
        Psat = _wilson_psat(const.Tcs[0], const.Pcs[0], const.omegas[0], T_K)
        is_gas = T_K >= const.Tcs[0] or pressure_pa <= Psat
        parent = flash_obj.gas if is_gas else flash_obj.liquid
        phase = parent.to(T=T_K, P=pressure_pa, zs=[1.0])
        phase.constants = parent.constants
        phase.correlations = parent.correlations
        rho = max(float(phase.rho_mass()), 0.001)
        mu = max(_phase_viscosity(phase, T_K, pressure_pa, const, corr), 1e-10)
        MW = float(const.MWs[0])
        sigma = _SIGMA_PLACEHOLDER if is_gas else _macleod_sugden_sigma(const, phase)
        return (1.0 if is_gas else 0.0, rho, rho, MW, MW, mu, mu, sigma)

    phases = list(res.phases)
    if len(phases) == 1:
        ph = phases[0]
        rho = max(float(ph.rho_mass()), 0.001)
        mu = max(_phase_viscosity(ph, T_K, pressure_pa, const, corr), 1e-10)
        MW = float(ph.MW())
        beta_g = 1.0 if _gas_like(ph, T_K, const) else 0.0
        sigma = _SIGMA_PLACEHOLDER if beta_g == 1.0 else _macleod_sugden_sigma(const, ph)
        return (beta_g, rho, rho, MW, MW, mu, mu, sigma)

    i_g = min(range(len(phases)), key=lambda i: phases[i].rho_mass())
    i_l = 1 - i_g if len(phases) == 2 else max(range(len(phases)), key=lambda i: phases[i].rho_mass())
    gas, liq = phases[i_g], phases[i_l]
    beta_g = float(res.betas[i_g])
    return (
        beta_g,
        max(float(gas.rho_mass()), 0.001),
        max(float(liq.rho_mass()), 0.001),
        float(gas.MW()),
        float(liq.MW()),
        max(_phase_viscosity(gas, T_K, pressure_pa, const, corr, True), 1e-10),
        max(_phase_viscosity(liq, T_K, pressure_pa, const, corr, False), 1e-10),
        _macleod_sugden_sigma(const, liq, gas),
    )


@lru_cache(maxsize=512)
def _hc_phase_compositions(
    component_names: tuple[str, ...],
    pressure_pa: float,
    temperature_c: float,
    zs: tuple[float, ...],
    eos_name: str = "PR",
) -> tuple[tuple[float, ...], tuple[float, ...], float]:
    """Return (gas_zs, liquid_zs, VF) for an HC-only PT flash.

    gas_zs and liquid_zs are mole fractions per HC component in the gas and
    liquid phases respectively.  When single-phase both tuples equal the feed.
    """
    flash_obj = _get_flash_obj(component_names, eos_name)
    T_K = temperature_c + 273.15
    try:
        res = _tp_flash(component_names, pressure_pa, temperature_c, zs, eos_name)
    except ZeroDivisionError:
        const  = flash_obj.constants
        Psat   = _wilson_psat(const.Tcs[0], const.Pcs[0], const.omegas[0], T_K)
        is_gas = T_K >= const.Tcs[0] or pressure_pa <= Psat
        return (zs, zs, 1.0 if is_gas else 0.0)

    VF = float(res.VF) if res.VF is not None else 0.0
    gas_ph = getattr(res, "gas", None)
    liq_ph = getattr(res, "liquid0", None)
    gas_zs    = tuple(float(y) for y in gas_ph.zs) if gas_ph is not None else zs
    liquid_zs = tuple(float(x) for x in liq_ph.zs) if liq_ph is not None else zs
    return (gas_zs, liquid_zs, VF)


@lru_cache(maxsize=2048)
def _flash_properties(
    component_names: tuple[str, ...],
    pressure_pa: float,
    temperature_c: float,
    zs: tuple[float, ...],
    eos_name: str = "PR",
) -> tuple[float, float, float, float, float, float]:
    """Return (rho kg/m³, mu Pa·s, Cp J/(kg·K), k W/(m·K), VF –, free_water_frac –).

    ``free_water_frac`` is the mole fraction of the feed that is liquid water
    (immiscible-water model).  It is 0.0 when no water component is present or
    when all water remains in the vapour phase.

    Free-water logic (immiscible-water model)
    -----------------------------------------
    When water is detected in ``component_names``:

    1. Water is removed from the component list and the HC mole fractions are
       re-normalised to sum to 1.
    2. The HC-only flash is performed with ``FlashVL`` at the same T and P.
    3. Water-phase split: the maximum amount of water vapour that the gas
       stream can carry is limited by the water saturation pressure
       (Wagner / IAPWS-IF97).  The remainder condenses as free liquid water.
    4. Bulk mixture density is the three-phase volumetric average
       (φ_gas·ρ_gas + φ_HCliq·ρ_HCliq + φ_water·ρ_water, PIPEPHASE-style);
       μ, Cp, k use the HC flash result.
    5. The overall vapour fraction (VF) includes water vapour in the gas phase.

    Single-component mixtures
    -------------------------
    ``thermo``'s Michelsen stability test divides by ``N-1 = 0`` for a single
    component.  When the ``ZeroDivisionError`` occurs the phase is determined
    via the Wilson K-value Psat estimate and the flash is retried with
    ``solution='g'`` or ``solution='l'``, which bypasses the stability test.

    Results are cached by (component_names, pressure_pa, temperature_c, zs, eos_name).
    Callers should round inputs before calling to maximise cache hit rates.
    """
    T_K = temperature_c + 273.15

    # ── free-water path ───────────────────────────────────────────────────────
    water_idx = find_water_index(component_names)
    if water_idx is not None:
        z_w = zs[water_idx]
        hc_names = tuple(n for i, n in enumerate(component_names) if i != water_idx)
        hc_zs_raw = tuple(z for i, z in enumerate(zs) if i != water_idx)
        sum_hc = sum(hc_zs_raw)

        if sum_hc < 1e-10 or not hc_names:
            # Pure water: return liquid water properties
            rho_wl = water_liquid_density_kg_m3(T_K)
            mu_wl  = water_liquid_viscosity_pa_s(T_K)
            return (rho_wl, mu_wl, 4182.0, 0.60, 0.0, z_w)

        hc_zs_norm = tuple(z / sum_hc for z in hc_zs_raw)

        # HC-only flash (recursive call; handles single-component HC too)
        hc_rho, hc_mu, hc_Cp, hc_k, VF_hc, _ = _flash_properties(
            hc_names, pressure_pa, temperature_c, hc_zs_norm, eos_name
        )

        # Water phase split (immiscible-water model)
        n_wl, n_wv = free_water_split(z_w, VF_hc, sum_hc, T_K, pressure_pa)

        # Three-phase density = total mass / total volume, per mole of feed.
        # Each HC phase uses its own molar mass (the vapour is much lighter
        # than the feed, the condensate much heavier); water vapour travels
        # in the gas at the gas molar volume.
        beta_g, rho_gas_hc, rho_liq_hc, MW_g, MW_l, _, _, _ = _hc_phase_split(
            hc_names, pressure_pa, temperature_c, hc_zs_norm, eos_name
        )
        rho_wl = water_liquid_density_kg_m3(T_K)

        n_gas_hc = beta_g * sum_hc
        n_liq_hc = (1.0 - beta_g) * sum_hc
        V_gas    = (n_gas_hc + n_wv) * MW_g / (1000.0 * rho_gas_hc)
        V_liq    = n_liq_hc * MW_l / (1000.0 * rho_liq_hc)
        V_wat    = n_wl * WATER_MW / (1000.0 * rho_wl)
        V_tot    = V_gas + V_liq + V_wat
        mass     = (n_gas_hc * MW_g + n_liq_hc * MW_l + (n_wl + n_wv) * WATER_MW) / 1000.0
        bulk_rho = mass / V_tot if V_tot > 1e-30 else hc_rho

        # Overall vapour fraction: HC vapour + water vapour (moles per mole of feed)
        bulk_VF = VF_hc * sum_hc + n_wv

        # Volume-fraction-weighted μ, Cp, k blending with free liquid water.
        # When n_wl ≈ 0 the correction vanishes; when WC is significant it
        # matters (water μ ≈ 1e-3 Pa·s >> gas μ ≈ 2e-5 Pa·s).
        if V_tot > 1e-30 and V_wat > 1e-30:
            mu_wl  = water_liquid_viscosity_pa_s(T_K)
            phi_hc = (V_gas + V_liq) / V_tot
            phi_w  = V_wat / V_tot
            bulk_mu = phi_hc * hc_mu + phi_w * mu_wl
            bulk_Cp = phi_hc * hc_Cp + phi_w * 4182.0
            bulk_k  = phi_hc * hc_k  + phi_w * 0.60
        else:
            bulk_mu, bulk_Cp, bulk_k = hc_mu, hc_Cp, hc_k

        return (
            max(bulk_rho, 0.001),
            max(bulk_mu, 1e-10),
            max(bulk_Cp, 1.0),
            max(bulk_k, 1e-6),
            bulk_VF,
            n_wl,
        )

    # ── standard HC-only path (no water) ─────────────────────────────────────
    flash_obj = _get_flash_obj(component_names, eos_name)
    try:
        res = _tp_flash(component_names, pressure_pa, temperature_c, zs, eos_name)
    except ZeroDivisionError:
        # thermo's Michelsen stability test divides by N-1=0 for single-component
        # mixtures.  Determine phase via Wilson Psat and evaluate the EOS phase
        # object directly — bypassing the flash entirely.
        const = flash_obj.constants
        Psat = _wilson_psat(const.Tcs[0], const.Pcs[0], const.omegas[0], T_K)
        parent = flash_obj.gas if (T_K >= const.Tcs[0] or pressure_pa <= Psat) else flash_obj.liquid
        VF_sc  = 1.0 if parent is flash_obj.gas else 0.0
        phase  = parent.to(T=T_K, P=pressure_pa, zs=[1.0])
        # .to() returns a new instance that does not inherit constants/correlations;
        # copy them from the cached parent so that rho_mass(), mu(), etc. work.
        phase.constants    = parent.constants
        phase.correlations = parent.correlations
        rho = max(float(phase.rho_mass()), 0.001)
        mu  = max(_phase_viscosity(phase, T_K, pressure_pa, const, parent.correlations), 1e-10)
        Cp  = max(float(phase.Cp_mass()),   1.0)
        k   = max(float(phase.k()),         1e-6)
        return (rho, mu, Cp, k, VF_sc, 0.0)
    except Exception as _exc:
        raise RuntimeError(
            f"EOS flash failed at T={temperature_c:.1f} °C, "
            f"P={pressure_pa / 1e6:.3f} MPa. "
            "The solver likely diverged to unphysical pressures. "
            "Check boundary conditions (pressure/flow BCs) and try reducing "
            "the relaxation factor in Numerics settings."
        ) from _exc

    rho = float(res.rho_mass())
    mu  = _result_viscosity(res, T_K, pressure_pa, flash_obj.constants, flash_obj.correlations)
    Cp  = float(res.Cp_mass())
    k   = float(res.k())
    VF  = float(res.VF) if res.VF is not None else 0.0

    return (max(rho, 0.001), max(mu, 1e-10), max(Cp, 1.0), max(k, 1e-6), VF, 0.0)


_CP_WATER = 4182.0


def _water_liquid_enthalpy(temperature_c: float, pressure_pa: float) -> float:
    """Liquid water, J/kg (0 °C, 1 atm reference): cp·T + (P − P_atm)/ρ."""
    return _CP_WATER * temperature_c + (pressure_pa - _ATM_PA) / water_liquid_density_kg_m3(temperature_c + 273.15)


def _water_latent_heat(temperature_c: float) -> float:
    """Latent heat of vaporisation of water, J/kg (linear fit, within ~1.5 %
    of steam tables from 0 to 150 °C)."""
    return 2.501e6 - 2369.0 * temperature_c


def _hc_enthalpy(component_names, pressure_pa, temperature_c, zs, eos_name) -> float:
    """Equilibrium specific enthalpy (J/kg) of an HC-only feed, from thermo."""
    flash_obj = _get_flash_obj(component_names, eos_name)
    T_K = temperature_c + 273.15
    try:
        res = _tp_flash(component_names, pressure_pa, temperature_c, zs, eos_name)
        return float(res.H()) / float(res.MW()) * 1000.0
    except ZeroDivisionError:
        const = flash_obj.constants
        Psat = _wilson_psat(const.Tcs[0], const.Pcs[0], const.omegas[0], T_K)
        parent = flash_obj.gas if (T_K >= const.Tcs[0] or pressure_pa <= Psat) else flash_obj.liquid
        phase = parent.to(T=T_K, P=pressure_pa, zs=[1.0])
        phase.constants = parent.constants
        phase.correlations = parent.correlations
        return float(phase.H()) / float(const.MWs[0]) * 1000.0


@lru_cache(maxsize=8192)
def _mixture_enthalpy(
    component_names: tuple[str, ...],
    pressure_pa: float,
    temperature_c: float,
    zs: tuple[float, ...],
    eos_name: str = "PR",
) -> float:
    """Equilibrium specific enthalpy of the feed, J/kg.

    Hydrocarbons: thermo's EOS flash enthalpy (includes the latent heat of
    condensation and the real-gas departure).  Free water (immiscible model):
    liquid water at cp·T + ΔP/ρ, water vapour carried by the gas at the liquid
    enthalpy plus the latent heat.  The reference state is arbitrary but fixed
    for a given composition, so enthalpy differences — all the energy balance
    uses — are meaningful.
    """
    water_idx = find_water_index(component_names)
    if water_idx is None:
        return _hc_enthalpy(component_names, pressure_pa, temperature_c, zs, eos_name)
    z_w = zs[water_idx]
    hc_names = tuple(n for i, n in enumerate(component_names) if i != water_idx)
    hc_zs_raw = tuple(z for i, z in enumerate(zs) if i != water_idx)
    sum_hc = sum(hc_zs_raw)
    h_wl = _water_liquid_enthalpy(temperature_c, pressure_pa)
    if sum_hc < 1e-10 or not hc_names:
        return h_wl
    hc_zs = tuple(z / sum_hc for z in hc_zs_raw)
    VF_hc = _flash_properties(hc_names, pressure_pa, temperature_c, hc_zs, eos_name)[4]
    n_wl, n_wv = free_water_split(z_w, VF_hc, sum_hc, temperature_c + 273.15, pressure_pa)
    MWs = _get_flash_obj(hc_names, eos_name).constants.MWs
    m_hc = sum_hc * sum(z * mw for z, mw in zip(hc_zs, MWs))
    m_wl, m_wv = n_wl * WATER_MW, n_wv * WATER_MW
    h_hc = _hc_enthalpy(hc_names, pressure_pa, temperature_c, hc_zs, eos_name)
    h_wv = h_wl + _water_latent_heat(temperature_c)
    return (m_hc * h_hc + m_wl * h_wl + m_wv * h_wv) / (m_hc + m_wl + m_wv)


@lru_cache(maxsize=1024)
def _reference_enthalpy(component_names: tuple[str, ...], zs: tuple[float, ...], eos_name: str = "PR") -> float:
    """Mixture enthalpy (J/kg) of the hydrocarbons as an ideal gas at 0 °C and
    the water as liquid at 0 °C (see ``_mixture_enthalpy`` for the basis)."""
    water_idx = find_water_index(component_names)
    if water_idx is not None:
        hc_names = tuple(n for i, n in enumerate(component_names) if i != water_idx)
        hc_raw = tuple(z for i, z in enumerate(zs) if i != water_idx)
        z_w = zs[water_idx]
    else:
        hc_names, hc_raw, z_w = component_names, zs, 0.0
    sum_hc = sum(hc_raw)
    if sum_hc < 1e-10 or not hc_names:
        return 0.0
    hc_zs = [z / sum_hc for z in hc_raw]
    flash_obj = _get_flash_obj(hc_names, eos_name)
    gas = flash_obj.gas.to(T=273.15, P=1.0, zs=hc_zs)          # 1 Pa: ideal gas
    gas.constants, gas.correlations = flash_obj.gas.constants, flash_obj.gas.correlations
    MW_hc = sum(z * mw for z, mw in zip(hc_zs, flash_obj.constants.MWs))
    h_hc = float(gas.H()) / MW_hc * 1000.0
    m_hc, m_w = sum_hc * MW_hc, z_w * WATER_MW
    return m_hc * h_hc / (m_hc + m_w)                            # water reference is 0


@lru_cache(maxsize=2048)
def _beggs_brill_properties(
    component_names: tuple[str, ...],
    pressure_pa: float,
    temperature_c: float,
    zs: tuple[float, ...],
    eos_name: str = "PR",
) -> tuple[float, float, float, float, float, float]:
    """Gas/liquid split for slip correlations: ``(C_L, rho_l, rho_g, mu_l, mu_g, sigma_l)``.

    The liquid is the HC condensate plus free water (volume-weighted density,
    viscosity and surface tension); the gas carries the water vapour.  C_L is
    the no-slip liquid volume fraction, and C_L·ρ_l + (1 − C_L)·ρ_g equals the
    bulk density of ``_flash_properties``.
    """
    T_K = temperature_c + 273.15
    water_idx = find_water_index(component_names)
    if water_idx is not None:
        z_w = zs[water_idx]
        hc_names = tuple(n for i, n in enumerate(component_names) if i != water_idx)
        hc_zs_raw = tuple(z for i, z in enumerate(zs) if i != water_idx)
        sum_hc = sum(hc_zs_raw)
        if sum_hc < 1e-10 or not hc_names:
            rho_w = water_liquid_density_kg_m3(T_K)
            return (1.0, rho_w, 1.2, water_liquid_viscosity_pa_s(T_K), 1.5e-5,
                    water_surface_tension(temperature_c, pressure_pa))
        hc_zs = tuple(z / sum_hc for z in hc_zs_raw)
        VF_hc = _flash_properties(hc_names, pressure_pa, temperature_c, hc_zs, eos_name)[4]
        n_wl, n_wv = free_water_split(z_w, VF_hc, sum_hc, T_K, pressure_pa)
    else:
        hc_names, hc_zs, sum_hc, n_wl, n_wv = component_names, zs, 1.0, 0.0, 0.0

    beta_g, rho_g, rho_l, MW_g, MW_l, mu_g, mu_l, sigma_l = _hc_phase_split(
        hc_names, pressure_pa, temperature_c, hc_zs, eos_name
    )
    n_g, n_l = beta_g * sum_hc, (1.0 - beta_g) * sum_hc
    V_gas = (n_g + n_wv) * MW_g / rho_g
    V_hcl = n_l * MW_l / rho_l
    rho_w = water_liquid_density_kg_m3(T_K)
    V_wat = n_wl * WATER_MW / rho_w
    V_liq = V_hcl + V_wat
    V_tot = V_gas + V_liq
    if V_tot <= 0.0:
        return (1.0, rho_l, rho_g, mu_l, mu_g, sigma_l)
    C_L = V_liq / V_tot
    if V_liq > 0.0:
        mu_w = water_liquid_viscosity_pa_s(T_K)
        sigma_w = water_surface_tension(temperature_c, pressure_pa)
        rho_liq = (V_hcl * rho_l + V_wat * rho_w) / V_liq
        mu_liq = (V_hcl * mu_l + V_wat * mu_w) / V_liq
        sig_liq = (V_hcl * sigma_l + V_wat * sigma_w) / V_liq
    else:
        rho_liq, mu_liq, sig_liq = rho_l, mu_l, sigma_l
    rho_gas = (n_g * MW_g + n_wv * WATER_MW) / V_gas if V_gas > 0.0 else rho_g
    return (C_L, rho_liq, rho_gas, mu_liq, mu_g, sig_liq)


class CompositionalFluid(FluidModel):
    """Equation-of-state fluid model backed by the ``thermo`` library.

    Performs a PT flash at each pipe's average pressure and temperature to
    compute density, viscosity, specific heat capacity, and thermal
    conductivity.  Two-phase mixtures are handled with the homogeneous
    no-slip model (volumetric-fraction-weighted mixture properties).

    **Free water**: if ``"water"`` (or ``"h2o"``, ``"7732-18-5"``) appears in
    ``components``, the solver uses the immiscible-water model automatically:
    the HC flash is performed on the water-free normalised composition, and
    free liquid water is quantified separately via the Wagner water saturation
    pressure.  No change to the API is required — include water in the
    component list and mole-fraction vector.  The phase-envelope is computed
    on a dry (water-free) basis as in HYSYS.

    The per-pipe mole-fraction vector is read from ``link_state.zs`` each
    time a property is evaluated.  ``SteadyCompositionalSolver`` propagates
    compositions from inlet boundary conditions through the network and writes
    the result to each ``PipeState.zs`` before calling the hydraulic solver.
    For link types that carry no ``zs`` attribute (fittings, pumps, heat
    sources), ``default_zs`` is used as a fallback.

    Args:
        components: Names of fluid components recognised by ``thermo``
            (e.g. ``["methane", "ethane", "water"]``).
        default_zs: Overall mole fractions used when a pipe's ``zs`` is not
            yet set.  Must sum to 1.0.
        eos_name: Equation of state key.  Supported values:

            - ``"PR"``     — Peng-Robinson (1976), default
            - ``"PR78"``   — Peng-Robinson (1978 alpha correction)
            - ``"PRSV"``   — Peng-Robinson-Stryjek-Vera (polar components)
            - ``"SRK"``    — Soave-Redlich-Kwong (1972)
            - ``"APISRK"`` — API-modified SRK (petroleum standard)
    """

    def __init__(
        self,
        components: list[str] | tuple[str, ...],
        default_zs: list[float] | tuple[float, ...],
        eos_name: str = "PR",
    ) -> None:
        self.component_names: tuple[str, ...] = tuple(components)
        self.default_zs: tuple[float, ...] = tuple(default_zs)
        eos_name = eos_name.upper()
        if eos_name not in SUPPORTED_EOS:
            raise ValueError(
                f"eos_name must be one of {list(SUPPORTED_EOS)} (got {eos_name!r})"
            )
        self.eos_name: str = eos_name
        if len(self.component_names) != len(self.default_zs):
            raise ValueError(
                f"components ({len(self.component_names)}) and default_zs "
                f"({len(self.default_zs)}) must have the same length"
            )
        if any(z < 0.0 for z in self.default_zs):
            raise ValueError("default_zs must be non-negative")
        s = sum(self.default_zs)
        if abs(s - 1.0) > 1e-6:
            raise ValueError(f"default_zs must sum to 1.0 (got {s:.6f})")

    # ── component molecular weights ───────────────────────────────────────────

    @cached_property
    def component_mws(self) -> tuple[float, ...]:
        """Molecular weights (g/mol) for each component, in the same order as component_names."""
        from thermo import ChemicalConstantsPackage
        constants, _ = ChemicalConstantsPackage.from_IDs(list(self.component_names))
        return tuple(constants.MWs)

    # ── helpers ───────────────────────────────────────────────────────────────

    def _get_zs(self, link_state) -> tuple[float, ...]:
        zs = getattr(link_state, "zs", None)
        return tuple(zs) if zs else self.default_zs

    @staticmethod
    def _get_pressure_pa(link_state) -> float:
        try:
            P_s = link_state.start_node.pressure_pa
            P_e = link_state.end_node.pressure_pa
            P_s = P_s if P_s is not None else 101325.0
            P_e = P_e if P_e is not None else 101325.0
            return max(0.5 * (P_s + P_e), 1.0)
        except AttributeError:
            return 101325.0

    @staticmethod
    def _get_temperature_c(link_state) -> float:
        t = getattr(link_state, "temperature_c", None)
        return float(t) if t is not None else 20.0

    def _props(self, link_state) -> tuple[float, float, float, float, float, float]:
        zs = self._get_zs(link_state)
        P  = self._get_pressure_pa(link_state)
        T  = self._get_temperature_c(link_state)
        P_r, T_r = _cache_key_pt(P, T)
        zs_r = tuple(round(z, _ZS_DIGITS) for z in zs)
        return _flash_properties(self.component_names, P_r, T_r, zs_r, self.eos_name)

    # ── FluidModel interface ──────────────────────────────────────────────────

    def density_for_link(self, link_state) -> float:
        return self._props(link_state)[0]

    def viscosity_for_link(self, link_state) -> float:
        return self._props(link_state)[1]

    def specific_heat_for_link(self, link_state) -> float:
        return self._props(link_state)[2]

    def thermal_conductivity_for_link(self, link_state) -> float:
        return self._props(link_state)[3]

    def vapor_fraction_for_link(self, link_state) -> float:
        """Return the equilibrium vapour fraction (VF) at the link's average P and T.

        Returns 0.0 for all-liquid, 1.0 for all-gas, and a value in (0, 1) for
        two-phase conditions.  When free water is present the VF includes both
        HC vapour and water vapour.  Reuses the cached flash result from _props.
        """
        try:
            return self._props(link_state)[4]
        except RuntimeError as exc:
            import warnings
            warnings.warn(str(exc), RuntimeWarning, stacklevel=2)
            return float("nan")

    def free_water_fraction_for_link(self, link_state) -> float:
        """Return the mole fraction of the feed that is liquid free water.

        Returns 0.0 when no water component is present or when all water
        remains in the vapour phase.  Values > 0 indicate that free liquid
        water is present at the link's average pressure and temperature.
        """
        try:
            return self._props(link_state)[5]
        except RuntimeError:
            return 0.0

    # ── enthalpy (energy balance: Joule-Thomson, latent heat) ─────────────────

    def _rounded(self, link_state, pressure_pa: float, temperature_c: float):
        P_r, T_r = _cache_key_pt(pressure_pa, temperature_c)
        zs_r = tuple(round(z, _ZS_DIGITS) for z in self._get_zs(link_state))
        return P_r, T_r, zs_r

    def enthalpy_for_link(self, link_state, pressure_pa: float, temperature_c: float) -> float:
        """Equilibrium specific enthalpy (J/kg) of the link's mixture at (P, T).

        Rounded far finer than the property cache (1 Pa, 1e-6 K): h enters the
        energy balance through differences, and a 1 mK step in T would leave a
        staircase that stalls the energy iteration."""
        zs_r = tuple(round(z, _ZS_DIGITS) for z in self._get_zs(link_state))
        P_r = float(round(pressure_pa))
        T_r = round(temperature_c, 6)
        return _mixture_enthalpy(self.component_names, P_r, T_r, zs_r, self.eos_name)

    def reference_enthalpy_for_link(self, link_state) -> float:
        """Enthalpy (J/kg, same basis as ``enthalpy_for_link``) of the link's
        mixture as an ideal gas at 0 °C — liquid water at 0 °C for free water —
        the reference of the global energy balance (additive by mass, so mixing
        streams does not shift it)."""
        zs_r = tuple(round(z, _ZS_DIGITS) for z in self._get_zs(link_state))
        return _reference_enthalpy(self.component_names, zs_r, self.eos_name)

    def specific_heat_at(self, link_state, pressure_pa: float, temperature_c: float) -> float:
        """Heat capacity (J/kg·K) of the link's mixture at (P, T), as in
        ``specific_heat_for_link``."""
        P_r, T_r, zs_r = self._rounded(link_state, pressure_pa, temperature_c)
        return _flash_properties(self.component_names, P_r, T_r, zs_r, self.eos_name)[2]

    def temperature_from_enthalpy_for_link(
        self, link_state, enthalpy_j_per_kg: float, pressure_pa: float, T_guess_c: float = 20.0,
    ) -> float:
        """Invert h(P, T) = enthalpy at constant composition (secant steps
        safeguarded by bisection, since h has kinks at phase boundaries)."""
        def f(T: float) -> float:
            return self.enthalpy_for_link(link_state, pressure_pa, T) - enthalpy_j_per_kg

        lo, hi = T_guess_c - 5.0, T_guess_c + 5.0
        f_lo, f_hi = f(lo), f(hi)
        for _ in range(40):                       # widen until the root is bracketed
            if f_lo * f_hi <= 0.0:
                break
            if abs(f_lo) < abs(f_hi):
                lo -= 2.0 * (hi - lo)
                f_lo = f(lo)
            else:
                hi += 2.0 * (hi - lo)
                f_hi = f(hi)
        else:
            return T_guess_c
        for _ in range(60):
            T = hi - f_hi * (hi - lo) / (f_hi - f_lo) if f_hi != f_lo else 0.5 * (lo + hi)
            if not (lo < T < hi):
                T = 0.5 * (lo + hi)
            f_T = f(T)
            if abs(f_T) < 1e-3 or hi - lo < 1e-5:
                return T
            if f_lo * f_T <= 0.0:
                hi, f_hi = T, f_T
                if abs(f_lo) > 2.0 * abs(f_T):        # Illinois step against stagnation
                    f_lo *= 0.5
            else:
                lo, f_lo = T, f_T
                if abs(f_hi) > 2.0 * abs(f_T):
                    f_hi *= 0.5
        return 0.5 * (lo + hi)

    def beggs_brill_phase_properties(self, link_state) -> tuple[float, float, float, float, float, float]:
        """``(C_L, rho_l, rho_g, mu_l, mu_g, sigma_l)`` at the link's mean P, T, z
        for slip correlations such as Beggs-Brill."""
        zs = self._get_zs(link_state)
        P_r, T_r = _cache_key_pt(self._get_pressure_pa(link_state), self._get_temperature_c(link_state))
        zs_r = tuple(round(z, _ZS_DIGITS) for z in zs)
        return _beggs_brill_properties(self.component_names, P_r, T_r, zs_r, self.eos_name)

    def flash_phase_detail_for_link(
        self, link_state
    ) -> tuple[float | None, tuple[float, ...], tuple[float, ...]]:
        """Return (liquid_fraction, gas_phase_zs, liquid_phase_zs) for the HC flash.

        liquid_fraction: mole fraction of feed that is HC liquid (feed basis).
        gas_phase_zs:    per-HC-component mole fractions in the gas phase.
        liquid_phase_zs: per-HC-component mole fractions in the HC liquid phase.

        Returns (None, (), ()) when no EOS flash is available.
        """
        try:
            zs = self._get_zs(link_state)
            P  = self._get_pressure_pa(link_state)
            T  = self._get_temperature_c(link_state)
            P_r, T_r = _cache_key_pt(P, T)

            water_idx = find_water_index(self.component_names)
            if water_idx is not None:
                hc_names   = tuple(n for i, n in enumerate(self.component_names) if i != water_idx)
                hc_zs_raw  = tuple(z for i, z in enumerate(zs) if i != water_idx)
                sum_hc     = sum(hc_zs_raw)
                if sum_hc < 1e-10 or not hc_names:
                    return (None, (), ())
                hc_zs_norm = tuple(round(z / sum_hc, _ZS_DIGITS) for z in hc_zs_raw)
                gas_zs, liquid_zs, VF_hc = _hc_phase_compositions(
                    hc_names, P_r, T_r, hc_zs_norm, self.eos_name
                )
                liq_frac = (1.0 - VF_hc) * sum_hc
            else:
                zs_r = tuple(round(z, _ZS_DIGITS) for z in zs)
                gas_zs, liquid_zs, VF_hc = _hc_phase_compositions(
                    self.component_names, P_r, T_r, zs_r, self.eos_name
                )
                liq_frac = 1.0 - VF_hc

            return (liq_frac, gas_zs, liquid_zs)
        except Exception:
            return (None, (), ())
