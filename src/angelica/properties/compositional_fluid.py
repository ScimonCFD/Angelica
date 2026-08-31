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

# Flash-object cache keyed by (component_names, eos_name) — avoids recreating
# ChemicalConstantsPackage + CEOSGas/Liquid on every cache miss.
_FLASH_OBJS: dict = {}


def _get_flash_obj(component_names: tuple[str, ...], eos_name: str):
    key = (component_names, eos_name)
    if key not in _FLASH_OBJS:
        import warnings

        from thermo import ChemicalConstantsPackage
        from thermo.eos_mix import PRMIX, SRKMIX
        from thermo.flash import FlashVL
        from thermo.phases import CEOSGas, CEOSLiquid

        from .phase_envelope import _build_kij_matrix
        eos_cls = SRKMIX if eos_name == "SRK" else PRMIX
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


def _wilson_psat(Tc: float, Pc: float, omega: float, T_K: float) -> float:
    """Wilson K-value estimate for saturation pressure (Pa) of a single component."""
    if T_K >= Tc:
        return Pc
    return Pc * math.exp(5.373 * (1.0 + omega) * (1.0 - Tc / T_K))


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
    4. Bulk mixture density accounts for the liquid-water volume contribution;
       all other bulk properties (μ, Cp, k) use the HC flash result.
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

        # Bulk density: volume-weighted mix of HC bulk and liquid water
        rho_wl = water_liquid_density_kg_m3(T_K)
        v_hc = sum_hc / hc_rho           # relative HC volume (arbitrary unit)
        v_wl = n_wl * WATER_MW / rho_wl  # relative liquid-water volume (same unit)
        phi_wl = v_wl / (v_hc + v_wl) if (v_hc + v_wl) > 1e-30 else 0.0
        bulk_rho = hc_rho * (1.0 - phi_wl) + rho_wl * phi_wl

        # Overall vapour fraction: HC vapour + water vapour (moles per mole of feed)
        bulk_VF = VF_hc * sum_hc + n_wv

        return (
            max(bulk_rho, 0.001),
            max(hc_mu, 1e-10),
            max(hc_Cp, 1.0),
            max(hc_k, 1e-6),
            bulk_VF,
            n_wl,
        )

    # ── standard HC-only path (no water) ─────────────────────────────────────
    flash_obj = _get_flash_obj(component_names, eos_name)
    try:
        res = flash_obj.flash(T=T_K, P=pressure_pa, zs=list(zs))
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
        mu  = max(float(phase.mu()),        1e-10)
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
    mu  = float(res.mu())
    Cp  = float(res.Cp_mass())
    k   = float(res.k())
    VF  = float(res.VF) if res.VF is not None else 0.0

    return (max(rho, 0.001), max(mu, 1e-10), max(Cp, 1.0), max(k, 1e-6), VF, 0.0)


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
        eos_name: Equation of state — ``"PR"`` (Peng-Robinson, default) or
            ``"SRK"`` (Soave-Redlich-Kwong).
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
        if eos_name not in ("PR", "SRK"):
            raise ValueError(f"eos_name must be 'PR' or 'SRK' (got {eos_name!r})")
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
        P_r  = float(round(P / 100.0) * 100.0)
        T_r  = round(T, 1)
        zs_r = tuple(round(z, 4) for z in zs)
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
