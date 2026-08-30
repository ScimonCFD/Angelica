from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

from .base import FluidModel


@dataclass(frozen=True)
class ThermalFluid(FluidModel):
    """Fluid model for non-isothermal simulation.

    All four properties can be supplied as either a constant (float) or a
    callable T_c -> value, where T_c is the local temperature in degrees
    Celsius read from link_state.temperature_c.  When a constant is given
    it is wrapped internally so the rest of the code sees a uniform interface.

    For the isothermal solver, use SingleComponentFluid instead.

    Args:
        reference_temperature_c: Fallback temperature (°C) used to evaluate
            fluid properties when a link state carries no temperature
            (e.g. Fittings and Pumps, which are hydraulic-only components).
            Defaults to the supply temperature of the first thermal BC if
            built via the solver; set explicitly when constructing manually.
    """

    _density_fn: Callable[[float], float]
    _viscosity_fn: Callable[[float], float]
    _specific_heat_fn: Callable[[float], float]
    _thermal_conductivity_fn: Callable[[float], float]
    reference_temperature_c: float = 20.0

    @staticmethod
    def _wrap(value: float | Callable[[float], float]) -> Callable[[float], float]:
        if callable(value):
            return value
        return lambda _t: float(value)

    @classmethod
    def from_constants(
        cls,
        density_kg_per_m3: float,
        viscosity_pa_s: float,
        specific_heat_j_per_kg_k: float,
        thermal_conductivity_w_per_m_k: float,
        reference_temperature_c: float = 20.0,
    ) -> ThermalFluid:
        return cls(
            _density_fn=cls._wrap(density_kg_per_m3),
            _viscosity_fn=cls._wrap(viscosity_pa_s),
            _specific_heat_fn=cls._wrap(specific_heat_j_per_kg_k),
            _thermal_conductivity_fn=cls._wrap(thermal_conductivity_w_per_m_k),
            reference_temperature_c=reference_temperature_c,
        )

    @classmethod
    def from_functions(
        cls,
        density_fn: Callable[[float], float],
        viscosity_fn: Callable[[float], float],
        specific_heat_fn: Callable[[float], float],
        thermal_conductivity_fn: Callable[[float], float],
        reference_temperature_c: float = 20.0,
    ) -> ThermalFluid:
        return cls(
            _density_fn=density_fn,
            _viscosity_fn=viscosity_fn,
            _specific_heat_fn=specific_heat_fn,
            _thermal_conductivity_fn=thermal_conductivity_fn,
            reference_temperature_c=reference_temperature_c,
        )

    def _temperature(self, link_state) -> float:
        t = getattr(link_state, "temperature_c", None)
        return float(t) if t is not None else self.reference_temperature_c

    def density_for_link(self, link_state) -> float:
        return self._density_fn(self._temperature(link_state))

    def viscosity_for_link(self, link_state) -> float:
        return self._viscosity_fn(self._temperature(link_state))

    def specific_heat_for_link(self, link_state) -> float:
        return self._specific_heat_fn(self._temperature(link_state))

    def thermal_conductivity_for_link(self, link_state) -> float:
        return self._thermal_conductivity_fn(self._temperature(link_state))
