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
from angelica.cases import (
    build_inline_heater_fixed_flow_case,
    build_symmetric_adiabatic_loop_case,
    build_symmetric_heat_loss_loop_case,
    build_thermal_mixing_junction_case,
)


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


def test_adiabatic_mixing_junction_matches_exact_mixing_temperature():
    """Two inlet streams should mix exactly at the junction in an adiabatic network."""
    result = SteadyNonIsothermalIncompressibleSolver().solve(
        build_thermal_mixing_junction_case()
    )

    assert result.converged
    expected_mix_c = 40.0
    assert result.node_temperatures_c[3] == pytest.approx(expected_mix_c, abs=0.1)
    assert result.node_temperatures_c[4] == pytest.approx(expected_mix_c, abs=0.1)

    flows = {cf.label: cf.mass_flow_kg_per_s for cf in result.component_flows}
    assert flows["Pipe:hot_branch"] == pytest.approx(1.0, abs=1e-6)
    assert flows["Pipe:cold_branch"] == pytest.approx(2.0, abs=1e-6)
    assert flows["Pipe:mixed_outlet"] == pytest.approx(3.0, abs=1e-6)


def test_inline_heater_fixed_flow_matches_exact_delta_t():
    """With fixed mass flow and adiabatic pipes, ΔT must match Q/(ṁcp)."""
    result = SteadyNonIsothermalIncompressibleSolver().solve(
        build_inline_heater_fixed_flow_case()
    )

    assert result.converged
    expected_delta_t = 50_000.0 / (1.0 * 4182.0)
    actual_delta_t = result.node_temperatures_c[4] - result.node_temperatures_c[1]
    assert actual_delta_t == pytest.approx(expected_delta_t, abs=0.2)

    flows = {cf.label: cf.mass_flow_kg_per_s for cf in result.component_flows}
    assert flows["Pipe:feed_pipe"] == pytest.approx(1.0, abs=1e-4)
    assert flows["HeatSource:heater"] == pytest.approx(1.0, abs=1e-4)
    assert flows["Pipe:exit_pipe"] == pytest.approx(1.0, abs=1e-4)


def test_symmetric_adiabatic_loop_stays_isothermal_and_splits_flow_evenly():
    result = SteadyNonIsothermalIncompressibleSolver().solve(
        build_symmetric_adiabatic_loop_case()
    )

    assert result.converged
    for node_id, temperature_c in result.node_temperatures_c.items():
        assert temperature_c == pytest.approx(60.0, abs=0.1), node_id

    flows = {cf.label: cf.mass_flow_kg_per_s for cf in result.component_flows}
    assert flows["Pipe:upper_branch"] == pytest.approx(1.0, abs=1e-4)
    assert flows["Pipe:lower_branch"] == pytest.approx(1.0, abs=1e-4)
    assert flows["Pipe:upper_return"] == pytest.approx(1.0, abs=1e-4)
    assert flows["Pipe:lower_return"] == pytest.approx(1.0, abs=1e-4)


def test_symmetric_heat_loss_loop_matches_exact_branch_temperature():
    result = SteadyNonIsothermalIncompressibleSolver().solve(
        build_symmetric_heat_loss_loop_case()
    )

    assert result.converged
    mdot_branch = 1.0
    cp = 4182.0
    U = 50.0
    D = 0.025
    L = 500.0
    T_in = 80.0
    T_amb = 20.0
    ntu = U * math.pi * D * L / (mdot_branch * cp)
    expected_branch_out = T_amb + (T_in - T_amb) * math.exp(-ntu)

    assert result.node_temperatures_c[3] == pytest.approx(expected_branch_out, abs=0.3)
    assert result.node_temperatures_c[4] == pytest.approx(expected_branch_out, abs=0.3)
    assert result.node_temperatures_c[5] == pytest.approx(expected_branch_out, abs=0.3)
    assert result.node_temperatures_c[6] == pytest.approx(expected_branch_out, abs=0.3)

    flows = {cf.label: cf.mass_flow_kg_per_s for cf in result.component_flows}
    assert flows["Pipe:upper_branch"] == pytest.approx(1.0, abs=1e-4)
    assert flows["Pipe:lower_branch"] == pytest.approx(1.0, abs=1e-4)
