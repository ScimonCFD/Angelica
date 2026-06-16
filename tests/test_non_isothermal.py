"""Tests for the non-isothermal incompressible solver.

Benchmark: single adiabatic-wall or heat-losing pipe against the analytical
exponential-decay solution for steady convection-diffusion with a uniform
heat-loss source term.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest

SRC_ROOT = Path(__file__).resolve().parents[1] / "src"
if str(SRC_ROOT) not in sys.path:
    sys.path.insert(0, str(SRC_ROOT))

from angelica import (
    NetworkCase,
    NonIsothermalSolverSettings,
    Pipe,
    PressureBoundary,
    SteadyNonIsothermalIncompressibleSolver,
    ThermalBoundary,
    ThermalFluid,
)
from angelica.closures import HybridScheme, PowerLawScheme, UpwindScheme


# ---------------------------------------------------------------------------
# Shared fixture: single pipe with heat loss
# ---------------------------------------------------------------------------
#
# Fluid: heavy-oil-like, constant properties
#   ρ = 900 kg/m³,  μ = 0.05 Pa·s,  cp = 2000 J/(kg·K),  k = 0.15 W/(m·K)
#
# Pipe: D = 0.01 m, L = 5 m, U = 100 W/(m²K), T_amb = 20 °C
#   Hagen-Poiseuille (Re ≈ 22, well laminar) gives
#     Q ≈ 9.82 × 10⁻⁶ m³/s  for ΔP = 10 000 Pa
#     ṁ ≈ 8.84 × 10⁻³ kg/s
#   NTU = U π D L / (ṁ cp) ≈ 0.89
#   T_out_analytical ≈ 44.7 °C  (inlet 80 °C, ambient 20 °C)
# ---------------------------------------------------------------------------

_FLUID = ThermalFluid.from_constants(
    density_kg_per_m3=900.0,
    viscosity_pa_s=0.05,
    specific_heat_j_per_kg_k=2000.0,
    thermal_conductivity_w_per_m_k=0.15,
)

_PIPE_WITH_HEAT_LOSS = Pipe(
    start_node=1,
    end_node=2,
    diameter_m=0.01,
    length_m=5.0,
    absolute_roughness_m=0.0,
    heat_transfer_coefficient_w_per_m2k=100.0,
    ambient_temperature_c=20.0,
    n_thermal_segments=50,
)

_PIPE_ADIABATIC = Pipe(
    start_node=1,
    end_node=2,
    diameter_m=0.01,
    length_m=5.0,
    absolute_roughness_m=0.0,
    heat_transfer_coefficient_w_per_m2k=0.0,
    ambient_temperature_c=20.0,
    n_thermal_segments=20,
)

_BOUNDARIES = dict(
    pressure_inlets=(PressureBoundary(node_id=1, pressure_pa=111_325.0),),
    pressure_outlets=(PressureBoundary(node_id=2, pressure_pa=101_325.0),),
    thermal_inlets=(ThermalBoundary(node_id=1, temperature_c=80.0),),
)


def _analytical_T_out(mdot: float) -> float:
    """Continuous-domain analytical outlet temperature (°C)."""
    U, D, L, cp = 100.0, 0.01, 5.0, 2000.0
    T_in, T_amb = 80.0, 20.0
    ntu = U * math.pi * D * L / (mdot * cp)
    return T_amb + (T_in - T_amb) * math.exp(-ntu)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_outlet_temperature_upwind():
    """Upwind scheme: outlet T within 1 °C of analytical exponential decay."""
    case = NetworkCase(
        name="heat_loss_upwind",
        fluid_model=_FLUID,
        components=(_PIPE_WITH_HEAT_LOSS,),
        **_BOUNDARIES,
    )
    solver = SteadyNonIsothermalIncompressibleSolver(
        non_isothermal_settings=NonIsothermalSolverSettings(
            max_temperature_iterations=30,
            temperature_tolerance_k=0.001,
        ),
        convection_scheme=UpwindScheme(),
    )
    result = solver.solve(case)

    assert result.converged

    mdot = result.component_flows[0].mass_flow_kg_per_s
    assert mdot > 0.0  # flow from node 1 → node 2

    T_out = result.node_temperatures_c[2]
    T_expected = _analytical_T_out(mdot)
    assert abs(T_out - T_expected) < 1.0, (
        f"outlet T = {T_out:.2f} °C, analytical = {T_expected:.2f} °C"
    )


@pytest.mark.parametrize("scheme", [HybridScheme(), PowerLawScheme()])
def test_outlet_temperature_alternative_schemes(scheme):
    """Hybrid and power-law schemes give similarly accurate outlet temperature."""
    case = NetworkCase(
        name="heat_loss_scheme",
        fluid_model=_FLUID,
        components=(_PIPE_WITH_HEAT_LOSS,),
        **_BOUNDARIES,
    )
    solver = SteadyNonIsothermalIncompressibleSolver(
        non_isothermal_settings=NonIsothermalSolverSettings(
            temperature_tolerance_k=0.001,
        ),
        convection_scheme=scheme,
    )
    result = solver.solve(case)

    assert result.converged
    mdot = result.component_flows[0].mass_flow_kg_per_s
    T_out = result.node_temperatures_c[2]
    T_expected = _analytical_T_out(mdot)
    assert abs(T_out - T_expected) < 1.5


def test_adiabatic_pipe_temperature_unchanged():
    """With U = 0 W/(m²K), outlet temperature equals inlet temperature."""
    case = NetworkCase(
        name="adiabatic",
        fluid_model=_FLUID,
        components=(_PIPE_ADIABATIC,),
        **_BOUNDARIES,
    )
    solver = SteadyNonIsothermalIncompressibleSolver()
    result = solver.solve(case)

    assert result.converged
    T_out = result.node_temperatures_c[2]
    assert abs(T_out - 80.0) < 0.5, f"adiabatic pipe: T_out = {T_out:.3f} °C, expected 80 °C"


def test_result_has_temperature_field():
    """SolveResult carries node_temperatures_c for all nodes."""
    case = NetworkCase(
        name="temp_field",
        fluid_model=_FLUID,
        components=(_PIPE_WITH_HEAT_LOSS,),
        **_BOUNDARIES,
    )
    result = SteadyNonIsothermalIncompressibleSolver().solve(case)
    assert set(result.node_temperatures_c.keys()) == {1, 2}
    assert result.node_temperatures_c[1] == pytest.approx(80.0, abs=0.1)


def test_no_flow_no_crash():
    """Zero ΔP (equal pressures) should converge without error."""
    case = NetworkCase(
        name="no_flow",
        fluid_model=_FLUID,
        components=(_PIPE_WITH_HEAT_LOSS,),
        pressure_inlets=(PressureBoundary(node_id=1, pressure_pa=101_325.0),),
        pressure_outlets=(PressureBoundary(node_id=2, pressure_pa=101_325.0),),
        thermal_inlets=(ThermalBoundary(node_id=1, temperature_c=80.0),),
    )
    result = SteadyNonIsothermalIncompressibleSolver().solve(case)
    # no assertion on T — just verify it doesn't raise
    assert result.node_temperatures_c is not None
