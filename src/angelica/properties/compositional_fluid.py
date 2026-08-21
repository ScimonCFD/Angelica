from __future__ import annotations

from functools import lru_cache

from .base import FluidModel


@lru_cache(maxsize=2048)
def _flash_properties(
    component_names: tuple[str, ...],
    pressure_pa: float,
    temperature_c: float,
    zs: tuple[float, ...],
) -> tuple[float, float, float, float]:
    """Return (rho kg/m³, mu Pa·s, Cp J/(kg·K), k W/(m·K)) for a mixture PT flash.

    Results are cached by (component_names, pressure_pa, temperature_c, zs).
    Callers should round inputs before calling to maximise cache hit rates.

    Two-phase handling uses the no-slip (homogeneous) model:
      - Volumetric fractions α_g, α_l from molar vapour fraction and phase MWs.
      - μ_mix = α_g·μ_g + α_l·μ_l
      - Cp_mix = (α_g·ρ_g·Cp_g + α_l·ρ_l·Cp_l) / ρ_mix
      - k_mix  = α_g·k_g  + α_l·k_l
    """
    from thermo import Mixture  # optional dependency

    T_K = temperature_c + 273.15
    try:
        m = Mixture(list(component_names), zs=list(zs), T=T_K, P=pressure_pa)
    except Exception as _exc:
        raise RuntimeError(
            f"EOS flash failed at T={temperature_c:.1f} °C, "
            f"P={pressure_pa / 1e6:.3f} MPa. "
            "The solver likely diverged to unphysical pressures. "
            "Check boundary conditions (pressure/flow BCs) and try reducing "
            "the relaxation factor in Numerics settings."
        ) from _exc
    rho = m.rho  # always available; homogeneous (no-slip) density

    if m.phase in ("g", "l", "s"):
        mu = m.mu
        Cp = m.Cp
        k  = m.k
    else:
        # Two-phase VLE: compute volumetric fractions from molar VF and phase MWs.
        VF   = m.VF
        L    = 1.0 - VF
        rhog = m.rhog or 1e-6
        rhol = m.rhol or 1e-6
        vol_g = VF * m.MWg / rhog
        vol_l = L  * m.MWl / rhol
        denom = vol_g + vol_l
        if denom < 1e-30:
            alpha_g, alpha_l = VF, L
        else:
            alpha_g = vol_g / denom
            alpha_l = 1.0 - alpha_g
        mu = alpha_g * (m.mug or 0.0) + alpha_l * (m.mul or 0.0)
        Cp = (
            (alpha_g * rhog * (m.Cpg or 0.0) + alpha_l * rhol * (m.Cpl or 0.0))
            / max(rho, 1e-30)
        )
        k  = alpha_g * (m.kg or 0.0) + alpha_l * (m.kl or 0.0)

    return (max(rho, 0.001), max(mu, 1e-10), max(Cp, 1.0), max(k, 1e-6))


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
    """

    def __init__(
        self,
        components: list[str] | tuple[str, ...],
        default_zs: list[float] | tuple[float, ...],
    ) -> None:
        self.component_names: tuple[str, ...] = tuple(components)
        self.default_zs: tuple[float, ...] = tuple(default_zs)
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

    # ── helpers ───────────────────────────────────────────────────────────────

    def _get_zs(self, link_state) -> tuple[float, ...]:
        zs = getattr(link_state, "zs", None)
        return tuple(zs) if zs else self.default_zs

    @staticmethod
    def _get_pressure_pa(link_state) -> float:
        try:
            P_s = link_state.start_node.pressure_pa or 101325.0
            P_e = link_state.end_node.pressure_pa or 101325.0
            return max(0.5 * (P_s + P_e), 1.0)
        except AttributeError:
            return 101325.0

    @staticmethod
    def _get_temperature_c(link_state) -> float:
        t = getattr(link_state, "temperature_c", None)
        return float(t) if t is not None else 20.0

    def _props(self, link_state) -> tuple[float, float, float, float]:
        zs = self._get_zs(link_state)
        P  = self._get_pressure_pa(link_state)
        T  = self._get_temperature_c(link_state)
        # Round inputs so nearby conditions share cache entries.
        P_r  = float(round(P / 100.0) * 100.0)
        T_r  = round(T, 1)
        zs_r = tuple(round(z, 4) for z in zs)
        return _flash_properties(self.component_names, P_r, T_r, zs_r)

    # ── FluidModel interface ──────────────────────────────────────────────────

    def density_for_link(self, link_state) -> float:
        return self._props(link_state)[0]

    def viscosity_for_link(self, link_state) -> float:
        return self._props(link_state)[1]

    def specific_heat_for_link(self, link_state) -> float:
        return self._props(link_state)[2]

    def thermal_conductivity_for_link(self, link_state) -> float:
        return self._props(link_state)[3]
