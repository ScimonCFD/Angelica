from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

from .base import FluidModel
from .eos import EquationOfState


@dataclass(frozen=True)
class CompressibleFluid(FluidModel):
    """Fluid model for compressible single-component gas simulation.

    Density is computed from an EquationOfState using the local pressure and
    temperature. Viscosity, specific heat, and thermal conductivity are
    treated as functions of temperature only (pressure correction is
    negligible for most gas pipeline applications).

    The local pressure is taken as the average of the link's inlet and outlet
    node pressures. If neither is available yet (early in the iteration),
    reference_pressure_pa is used as a fallback.

    Args:
        eos: Equation of state that provides density(pressure_pa, temperature_c).
        reference_pressure_pa: Fallback pressure (Pa) when node pressures are
            not yet initialised. Defaults to 101 325 Pa (1 atm).
        reference_temperature_c: Fallback temperature (°C) when no temperature
            is stored on the link state. Defaults to 20 °C.
    """

    eos: EquationOfState
    _viscosity_fn: Callable[[float], float]
    _specific_heat_fn: Callable[[float], float]
    _thermal_conductivity_fn: Callable[[float], float]
    reference_pressure_pa: float = 101_325.0
    reference_temperature_c: float = 20.0

    @staticmethod
    def _wrap(value: float | Callable[[float], float]) -> Callable[[float], float]:
        if callable(value):
            return value
        return lambda _t: float(value)

    @classmethod
    def from_constants(
        cls,
        eos: EquationOfState,
        viscosity_pa_s: float,
        specific_heat_j_per_kg_k: float,
        thermal_conductivity_w_per_m_k: float,
        reference_pressure_pa: float = 101_325.0,
        reference_temperature_c: float = 20.0,
    ) -> CompressibleFluid:
        return cls(
            eos=eos,
            _viscosity_fn=cls._wrap(viscosity_pa_s),
            _specific_heat_fn=cls._wrap(specific_heat_j_per_kg_k),
            _thermal_conductivity_fn=cls._wrap(thermal_conductivity_w_per_m_k),
            reference_pressure_pa=reference_pressure_pa,
            reference_temperature_c=reference_temperature_c,
        )

    @classmethod
    def from_functions(
        cls,
        eos: EquationOfState,
        viscosity_fn: Callable[[float], float],
        specific_heat_fn: Callable[[float], float],
        thermal_conductivity_fn: Callable[[float], float],
        reference_pressure_pa: float = 101_325.0,
        reference_temperature_c: float = 20.0,
    ) -> CompressibleFluid:
        return cls(
            eos=eos,
            _viscosity_fn=viscosity_fn,
            _specific_heat_fn=specific_heat_fn,
            _thermal_conductivity_fn=thermal_conductivity_fn,
            reference_pressure_pa=reference_pressure_pa,
            reference_temperature_c=reference_temperature_c,
        )

    def _pressure(self, link_state) -> float:
        p_start = getattr(link_state.start_node, "pressure_pa", None)
        p_end = getattr(link_state.end_node, "pressure_pa", None)
        if p_start is not None and p_end is not None:
            return 0.5 * (p_start + p_end)
        return p_start if p_start is not None else (p_end if p_end is not None else self.reference_pressure_pa)

    def _temperature(self, link_state) -> float:
        t = getattr(link_state, "temperature_c", None)
        return float(t) if t is not None else self.reference_temperature_c

    def density_for_link(self, link_state) -> float:
        return self.eos.density(self._pressure(link_state), self._temperature(link_state))

    def viscosity_for_link(self, link_state) -> float:
        return self._viscosity_fn(self._temperature(link_state))

    def specific_heat_for_link(self, link_state) -> float:
        return self._specific_heat_fn(self._temperature(link_state))

    def thermal_conductivity_for_link(self, link_state) -> float:
        return self._thermal_conductivity_fn(self._temperature(link_state))
