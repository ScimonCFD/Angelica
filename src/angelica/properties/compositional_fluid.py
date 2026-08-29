from __future__ import annotations

from functools import cached_property, lru_cache

from .base import FluidModel

# Flash-object cache keyed by (component_names, eos_name) — avoids recreating
# ChemicalConstantsPackage + CEOSGas/Liquid on every cache miss.
_FLASH_OBJS: dict = {}


def _get_flash_obj(component_names: tuple[str, ...], eos_name: str):
    key = (component_names, eos_name)
    if key not in _FLASH_OBJS:
        import warnings
        from thermo import ChemicalConstantsPackage
        from thermo.flash import FlashVL
        from thermo.phases import CEOSGas, CEOSLiquid
        from thermo.eos_mix import PRMIX, SRKMIX
        from .phase_envelope import _build_kij_matrix
        eos_cls = SRKMIX if eos_name == "SRK" else PRMIX
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            constants, props = ChemicalConstantsPackage.from_IDs(list(component_names))
            kijs = _build_kij_matrix(constants)
            eos_kw = dict(Tcs=constants.Tcs, Pcs=constants.Pcs, omegas=constants.omegas, kijs=kijs)
            gas_phase = CEOSGas(eos_cls, eos_kw, HeatCapacityGases=props.HeatCapacityGases)
            liq_phase = CEOSLiquid(eos_cls, eos_kw, HeatCapacityGases=props.HeatCapacityGases)
            flash_obj = FlashVL(constants, props, liquid=liq_phase, gas=gas_phase)
        _FLASH_OBJS[key] = flash_obj
    return _FLASH_OBJS[key]


@lru_cache(maxsize=2048)
def _flash_properties(
    component_names: tuple[str, ...],
    pressure_pa: float,
    temperature_c: float,
    zs: tuple[float, ...],
    eos_name: str = "PR",
) -> tuple[float, float, float, float, float]:
    """Return (rho kg/m³, mu Pa·s, Cp J/(kg·K), k W/(m·K), VF –) for a mixture PT flash.

    Results are cached by (component_names, pressure_pa, temperature_c, zs, eos_name).
    Callers should round inputs before calling to maximise cache hit rates.

    Two-phase handling: EquilibriumState bulk properties use volumetric-fraction
    mixing automatically via the thermo library.
    """
    flash_obj = _get_flash_obj(component_names, eos_name)
    T_K = temperature_c + 273.15
    try:
        res = flash_obj.flash(T=T_K, P=pressure_pa, zs=list(zs))
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

    return (max(rho, 0.001), max(mu, 1e-10), max(Cp, 1.0), max(k, 1e-6), VF)


class CompositionalFluid(FluidModel):
    """Equation-of-state fluid model backed by the ``thermo`` library.

    Performs a PT flash at each pipe's average pressure and temperature to
    compute density, viscosity, specific heat capacity, and thermal
    conductivity.  Two-phase mixtures are handled with the homogeneous
    no-slip model.

    The per-pipe mole-fraction vector is read from ``link_state.zs`` each
    time a property is evaluated.  ``SteadyCompositionalSolver`` propagates
    compositions from inlet boundary conditions through the network and writes
    the result to each ``PipeState.zs`` before calling the hydraulic solver.
    For link types that carry no ``zs`` attribute (fittings, pumps, heat
    sources), ``default_zs`` is used as a fallback.

    Args:
        components: Names of fluid components recognised by ``thermo``
            (e.g. ``["methane", "ethane", "propane"]``).
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

    def _props(self, link_state) -> tuple[float, float, float, float, float]:
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
        """Return the equilibrium vapor fraction (VF) at the link's average P and T.

        Returns 0.0 for all-liquid, 1.0 for all-gas, and a value in (0, 1) for
        two-phase conditions.  Reuses the cached flash result from _props.
        """
        try:
            return self._props(link_state)[4]
        except RuntimeError as exc:
            import warnings
            warnings.warn(str(exc), RuntimeWarning, stacklevel=2)
            return float("nan")
